#!/usr/bin/env bash
#SBATCH --job-name=catalog-qc
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/catalog-qc-%j.log
set -euo pipefail
module reset
export XDG_CACHE_HOME="$PWD/.cache"
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:/vast/palmer/apps/avx2/software/miniconda/24.7.1/condabin:$PATH"
record="data/results/catalog_qc/2026-10-20/executions/$SLURM_JOB_ID"
mkdir -p "$record"
exec > >(tee "$record/execution.log") 2>&1
trap 'code=$?; date -u > "$record/finished.txt"; printf "%s\n" "$code" > "$record/exit_code.txt"' EXIT
date -u > "$record/started.txt"
git rev-parse HEAD > "$record/git_head.txt"
git status --short > "$record/git_status.txt"
tar -czf "$record/source.tar.gz" config workflow scripts profiles envs
sha256sum "$record/source.tar.gz" > "$record/source_sha256.txt"
set -x
python -m snakemake --snakefile workflow/catalog_qc.smk domain_screen \
  --profile profiles/slurm --jobs 4 --nolock --rerun-triggers mtime \
  --default-resources partition=ycga runtime=120 mem_mb=4000
