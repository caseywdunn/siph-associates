#!/usr/bin/env bash
#SBATCH --job-name=deposit-validation-tools
#SBATCH --partition=ycga
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/deposit-validation-tools-%j.log
set -euo pipefail
module reset
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/barrnap_env/bin:$PATH"
repo="$PWD"
tools_root="$repo/data/software/deposit_validation"
record="$repo/data/results/deposit/executions/validation-tools-$SLURM_JOB_ID"
export TMPDIR="$repo/tmp/validation-tools-$SLURM_JOB_ID"
mkdir -p "$tools_root" "$record" "$TMPDIR"
cp scripts/install_deposit_validation_tools.sbatch "$record/"
trap 'status=$?; printf "%s\n" "$status" > "$record/exit_code.txt"' EXIT
# Isolated tool directories; the project environments and user shell stay unchanged.
if [[ ! -x "$tools_root/table2asn/table2asn" ]]; then
  mkdir -p "$tools_root/table2asn"
  curl --fail --location --retry 3 https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz \
    --output "$tools_root/table2asn/linux64.table2asn.gz.part"
  mv "$tools_root/table2asn/linux64.table2asn.gz.part" "$tools_root/table2asn/linux64.table2asn.gz"
  gzip -dc "$tools_root/table2asn/linux64.table2asn.gz" > "$tools_root/table2asn/table2asn.part"
  chmod u+x "$tools_root/table2asn/table2asn.part"
  mv "$tools_root/table2asn/table2asn.part" "$tools_root/table2asn/table2asn"
fi
"$tools_root/table2asn/table2asn" -help > "$record/table2asn_help.txt" 2>&1
sha256sum "$tools_root/table2asn/table2asn" > "$record/table2asn.sha256"
if [[ ! -x "$tools_root/ribovore/ribovore/ribotyper" || ! -x "$tools_root/ribovore/infernal/binaries/cmsearch" ]]; then
  mkdir -p "$tools_root/ribovore"
  cd "$tools_root/ribovore"
  curl --fail --location --retry 3 https://raw.githubusercontent.com/ncbi/ribovore/master/install.sh \
    --output install.sh.part
  mv install.sh.part install.sh
  # Review again if the upstream installer changes the release to be installed.
  grep -q 'VERSION="1.0.5"' install.sh
  cp install.sh "$record/ribovore_install.sh"
  sha256sum install.sh > "$record/ribovore_installer.sha256"
  bash install.sh rt-linux keep > "$record/ribovore_install.log" 2>&1
  test -x ribovore/ribotyper
  test -x infernal/binaries/cmsearch
fi
export RIBOINSTALLDIR="$tools_root/ribovore"
export RIBOSCRIPTSDIR="$RIBOINSTALLDIR/ribovore"
export RIBOINFERNALDIR="$RIBOINSTALLDIR/infernal/binaries"
export RIBOEASELDIR="$RIBOINFERNALDIR"
export RIBOSEQUIPDIR="$RIBOINSTALLDIR/sequip"
export RIBOTIMEDIR=/usr/bin
export PERL5LIB="$RIBOSCRIPTSDIR:$RIBOSEQUIPDIR:${PERL5LIB:-}"
perl -MTime::HiRes -e 'print "$^V\n"' > "$record/perl_version.txt"
"$RIBOSCRIPTSDIR/ribotyper" -h > "$record/ribotyper_help.txt" 2>&1
sha256sum "$RIBOSCRIPTSDIR/ribotyper" "$RIBOINFERNALDIR/cmsearch" > "$record/ribotyper.sha256"
"$RIBOSCRIPTSDIR/ribotyper" "$RIBOSCRIPTSDIR/testfiles/example-16.fa" "$record/ribotyper_example" \
  > "$record/ribotyper_example.log" 2>&1

touch "$tools_root/ribovore/installed.ok"
