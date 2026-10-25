#!/usr/bin/env bash
#SBATCH --job-name=deposit-checks
#SBATCH --partition=ycga
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/deposit-checks-%j.log
set -euo pipefail
module reset
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:$PATH"
export TMPDIR="$PWD/data/results/deposit/tmp/checks-$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
python -m pytest -q tests/test_deposit_preparation.py
python -m snakemake -s workflow/deposit_preparation.smk deposit_preparation --cores 16 --resources mem_mb=64000 --dry-run
singularity exec /gpfs/ycga/work/dunn/cwd7/databases/fcs/fcs-gx.sif /app/bin/run_gx --help
