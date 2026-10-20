#!/usr/bin/env bash
#SBATCH --job-name=assembly-completion
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/assembly-completion-%j.log
set -euo pipefail
module reset
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:/vast/palmer/apps/avx2/software/miniconda/24.7.1/condabin:$PATH"
export XDG_CACHE_HOME="$PWD/.cache"
record="data/results/assembly_completion/2026-10-20/executions/$SLURM_JOB_ID"
mkdir -p "$record"
exec > >(tee "$record/execution.log") 2>&1
trap 'code=$?; date -u > "$record/finished.txt"; printf "%s\n" "$code" > "$record/exit_code.txt"' EXIT
date -u > "$record/started.txt"
git rev-parse HEAD > "$record/git_head.txt"
git status --short > "$record/git_status.txt"
tar -czf "$record/source.tar.gz" config workflow scripts profiles envs
sha256sum "$record/source.tar.gz" > "$record/source_sha256.txt"
set -x
python -m snakemake --snakefile workflow/assembly_completion.smk assemblies \
  --profile profiles/slurm --jobs 6 --cores 10000 --local-cores 1 --nolock \
  --default-resources partition=ycga runtime=120 mem_mb=2000 \
  --rerun-triggers mtime --rerun-incomplete
