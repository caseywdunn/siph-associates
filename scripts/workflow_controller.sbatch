#!/usr/bin/env bash
#SBATCH --job-name=sa.controller
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm/controller.%j.out
set -euo pipefail

PROJECT="${SIPH_ASSOCIATES_PROJECT:?missing SIPH_ASSOCIATES_PROJECT}"
TARGET="${SIPH_ASSOCIATES_TARGET:?missing SIPH_ASSOCIATES_TARGET}"
CONFIG="${SIPH_ASSOCIATES_CONFIG:?missing SIPH_ASSOCIATES_CONFIG}"
SNAKEMAKE=/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/snakemake

cd "$PROJECT"
export XDG_CACHE_HOME="$PROJECT/.cache"
mkdir -p "$XDG_CACHE_HOME"
exec "$SNAKEMAKE" \
    "$TARGET" \
    --profile profiles/slurm \
    --configfile "$CONFIG" \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 180 \
    --rerun-triggers mtime
