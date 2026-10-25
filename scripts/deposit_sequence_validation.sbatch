#!/usr/bin/env bash
#SBATCH --job-name=deposit-sequence-validation
#SBATCH --partition=ycga
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/deposit-sequence-validation-%j.log
set -euo pipefail
module reset
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/barrnap_env/bin:$PATH"
python_bin=/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python
repo="$PWD"
record="$repo/data/results/deposit/executions/sequence-validation-$SLURM_JOB_ID"
export TMPDIR="$repo/tmp/sequence-validation-$SLURM_JOB_ID"
mkdir -p "$record" "$TMPDIR"
trap 'status=$?; printf "%s\n" "$status" > "$record/exit_code.txt"' EXIT
cp scripts/deposit_sequence_validation.sbatch workflow/scripts/run_table2asn_draft_checks.py workflow/scripts/package_deposit_sequences.py workflow/scripts/prepare_deposit_mag.py "$record/"
sha256sum data/results/deposit/bundle/sets/ssu.fna data/results/cnidarian_16s/mycoplasmatales_mag_16s.fasta \
  data/results/deposit/bundle/genbank_candidates.tsv > "$record/input.sha256"
export RIBOINSTALLDIR="$repo/data/software/deposit_validation/ribovore"
export RIBOSCRIPTSDIR="$RIBOINSTALLDIR/ribovore"
export RIBOINFERNALDIR="$RIBOINSTALLDIR/infernal/binaries"
export RIBOEASELDIR="$RIBOINFERNALDIR"
export RIBOSEQUIPDIR="$RIBOINSTALLDIR/sequip"
export RIBOTIMEDIR=/usr/bin
export PERL5LIB="$RIBOSCRIPTSDIR:$RIBOSEQUIPDIR:${PERL5LIB:-}"
"$RIBOSCRIPTSDIR/ribotyper" data/results/deposit/bundle/sets/ssu.fna "$record/eukaryote_ssu" \
  > "$record/eukaryote_ssu.log" 2>&1
"$RIBOSCRIPTSDIR/ribotyper" data/results/cnidarian_16s/mycoplasmatales_mag_16s.fasta "$record/mag_16s" \
  > "$record/mag_16s.log" 2>&1
"$python_bin" workflow/scripts/run_table2asn_draft_checks.py --bundle data/results/deposit/bundle \
  --table2asn "$repo/data/software/deposit_validation/table2asn/table2asn" \
  --output "$record/table2asn_drafts"
