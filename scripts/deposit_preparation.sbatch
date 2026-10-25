#!/usr/bin/env bash
#SBATCH --job-name=deposit-preparation
#SBATCH --partition=ycga
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/deposit-preparation-%j.log
set -euo pipefail
module reset
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:$PATH"
export TMPDIR="$PWD/data/results/deposit/tmp/controller-$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
record="$PWD/data/results/deposit/executions/preparation-$SLURM_JOB_ID"
mkdir -p "$record"
git rev-parse HEAD > "$record/git_head.txt"
git diff > "$record/source.diff"
tar -czf "$record/source.tar.gz" config/deposit_preparation.json config/mag_*.tsv workflow/deposit_preparation.smk workflow/rules/deposit_preparation.smk workflow/scripts/prepare_deposit_mag.py workflow/scripts/package_deposit_sequences.py workflow/scripts/prepare_deposit_validation_inputs.py config/viral_catalog_classifications.tsv figures/manuscript_revision/catalog/viral_associate.manifest.tsv figures/manuscript_revision/catalog/viral_endogenous_candidate.manifest.tsv
sha256sum "$record/source.tar.gz" > "$record/source.sha256"
trap 'status=$?; printf "%s\n" "$status" > "$record/exit_code.txt"' EXIT
python -m snakemake -s workflow/deposit_preparation.smk deposit_preparation \
  --cores 16 --resources mem_mb=64000 --rerun-triggers mtime --rerun-incomplete --printshellcmds "$@"
