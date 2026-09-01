#!/usr/bin/env bash
#SBATCH --job-name=sa.controller
#SBATCH --partition=week
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/slurm/controller.%j.out
set -euo pipefail

PROJECT="${SIPH_ASSOCIATES_PROJECT:?missing SIPH_ASSOCIATES_PROJECT}"
TARGET="${SIPH_ASSOCIATES_TARGET:?missing SIPH_ASSOCIATES_TARGET}"
CONFIG="${SIPH_ASSOCIATES_CONFIG:?missing SIPH_ASSOCIATES_CONFIG}"
SNAKEMAKE=/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/snakemake
CONDA_BIN=/vast/palmer/apps/avx2/software/miniconda/24.7.1/condabin

cd "$PROJECT"
export XDG_CACHE_HOME="$PROJECT/.cache"
export CONDA_PKGS_DIRS="$PROJECT/.cache/conda/pkgs"
export PATH="$CONDA_BIN:$PATH"
mkdir -p "$XDG_CACHE_HOME" "$CONDA_PKGS_DIRS"
# Full-cohort failures reached the original 24 GB allocation exactly. Keep
# these scheduler-only recovery overrides outside config/config.yaml so a
# resource correction does not invalidate accepted analytical outputs through
# the frozen run snapshot.
exec "$SNAKEMAKE" \
    "$TARGET" \
    --profile profiles/slurm \
    --configfile "$CONFIG" \
    --set-resources \
        trim_reads:mem_mb=48000 \
        phyloflash_screen:mem_mb=48000 \
        phyloflash_screen:runtime=360 \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 180 \
    --rerun-triggers mtime
