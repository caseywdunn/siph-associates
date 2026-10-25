#!/usr/bin/env bash
#SBATCH --job-name=revision-16s
#SBATCH --partition=ycga
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/revision-16s-%j.log
set -euo pipefail
module reset
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:$PATH"
export TMPDIR="$PWD/data/results/deposit/tmp/16s-$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
record="$PWD/data/results/deposit/executions/16s-$SLURM_JOB_ID"
mkdir -p "$record/baseline"
baseline="$PWD/data/results/deposit/baselines/16s-2026-10-24"
if [[ ! -d "$baseline" ]]; then
  mkdir -p "$baseline.part"
  for file in mycoplasmatales_mag_16s.fasta asv_matches.tsv matching_asv_samples.tsv matching_asv_hosts.tsv; do
    cp "data/results/cnidarian_16s/$file" "$baseline.part/$file"
  done
  mv "$baseline.part" "$baseline"
fi
cp "$baseline/"* "$record/baseline/"
sha256sum "$record/baseline/"* > "$record/baseline.sha256"
git rev-parse HEAD > "$record/git_head.txt"
git diff > "$record/source.diff"
cp config/mag_contig_removals.tsv workflow/scripts/extract_mag_16s.py workflow/scripts/validate_16s_revision.py workflow/scripts/prepare_deposit_mag.py workflow/rules/cnidarian_16s.smk "$record/"
find "$record" -maxdepth 1 -type f ! -name source.sha256 -print0 | sort -z | xargs -0 sha256sum > "$record/source.sha256"
trap 'status=$?; printf "%s\n" "$status" > "$record/exit_code.txt"' EXIT
python -m snakemake -s workflow/cnidarian_16s_revision.smk revision_16s --cores 8 \
  --resources mem_mb=64000 --rerun-triggers mtime --rerun-incomplete --printshellcmds "$@"

if [[ " $* " != *" --dry-run "* && " $* " != *" -n "* ]]; then
  python workflow/scripts/validate_16s_revision.py --before "$record/baseline" \
    --after data/results/cnidarian_16s --report "$record/before_after.json"
  cp data/results/cnidarian_16s/matching_asv_samples.tsv tables/mycoplasmatales/cnidarian_16s_matching_asv_samples.tsv
  cp data/results/cnidarian_16s/matching_asv_hosts.tsv tables/mycoplasmatales/cnidarian_16s_matching_asv_hosts.tsv
fi
