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
# A space-separated TARGET requests several targets in one DAG.
read -ra TARGETS <<< "$TARGET"
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
# Cohort aggregation likewise reached its original 2 GB allocation while
# retaining nomination rows from all 205 libraries; give the finishing job
# sufficient headroom without changing any analytical input or parameter.
exec "$SNAKEMAKE" \
    "${TARGETS[@]}" \
    --profile profiles/slurm \
    --configfile "$CONFIG" \
    --set-resources \
        trim_reads:mem_mb=48000 \
        phyloflash_screen:mem_mb=48000 \
        phyloflash_screen:runtime=360 \
        aggregate_phase2:mem_mb=8000 \
        screen_cohort:mem_mb=8000 \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 180 \
    --rerun-triggers mtime
