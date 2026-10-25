#!/usr/bin/env bash
#SBATCH --job-name=deposit-fcs-preflight
#SBATCH --partition=ycga
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/deposit-fcs-preflight-%j.log
set -euo pipefail
module reset
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:$PATH"
export TMPDIR="$PWD/tmp/fcs-preflight-$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
python -m snakemake -s workflow/deposit_preparation.smk check_deposit_fcsgx_inputs \
  --cores 8 --resources mem_mb=32000 --rerun-triggers mtime --rerun-incomplete --printshellcmds
python -m snakemake -s workflow/deposit_preparation.smk deposit_fcsgx \
  --cores 16 --resources mem_mb=720000 --config enable_fcsgx=true --dry-run
singularity exec /gpfs/ycga/work/dunn/cwd7/databases/fcs/fcs-gx.sif /app/bin/run_gx --help
