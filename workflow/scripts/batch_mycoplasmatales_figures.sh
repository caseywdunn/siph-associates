#!/bin/bash
#SBATCH --job-name=mycoplasmatales-figures
#SBATCH --partition=ycga
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=figures/mycoplasmatales/slurm-%j.log

set -euo pipefail
module reset
export MPLCONFIGDIR="${SLURM_SUBMIT_DIR}/data/results/manuscript_revision/2026-10-24/tmp/mycoplasmatales-mpl-${SLURM_JOB_ID}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
mkdir -p "$MPLCONFIGDIR"
trap 'rm -rf "$MPLCONFIGDIR"' EXIT
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python \
  workflow/scripts/plot_mycoplasmatales_figures.py
