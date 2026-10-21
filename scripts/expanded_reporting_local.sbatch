#!/usr/bin/env bash
#SBATCH --job-name=expanded-local-reporting
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/expanded-local-reporting-%j.log
set -euo pipefail
module reset
export XDG_CACHE_HOME="$PWD/.cache"
export MPLCONFIGDIR=/tmp/siph-associates-matplotlib
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:/vast/palmer/apps/avx2/software/miniconda/24.7.1/condabin:$PATH"
record="data/results/post_expansion/2026-10-21/reporting/executions/$SLURM_JOB_ID"
mkdir -p "$record"
exec > >(tee "$record/execution.log") 2>&1
trap 'code=$?; date -u > "$record/finished.txt"; printf "%s\n" "$code" > "$record/exit_code.txt"' EXIT
date -u > "$record/started.txt"
git rev-parse HEAD > "$record/git_head.txt"
git status --short > "$record/git_status.txt"
tar -czf "$record/source.tar.gz" config workflow scripts profiles envs figures/*.py
sha256sum "$record/source.tar.gz" > "$record/source_sha256.txt"
python -m snakemake -s workflow/expanded_reporting.smk "$@" \
  --cores "$SLURM_CPUS_PER_TASK" --resources mem_mb=128000 --use-conda --conda-frontend mamba \
  --nolock --rerun-incomplete --rerun-triggers mtime --printshellcmds
