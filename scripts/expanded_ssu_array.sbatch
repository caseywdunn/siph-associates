#!/usr/bin/env bash
#SBATCH --job-name=expanded-ssu-batch
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --array=1-6%6
#SBATCH --output=logs/expanded-ssu-%A_%a.log
set -euo pipefail
module reset
export XDG_CACHE_HOME="$PWD/.cache"
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:/vast/palmer/apps/avx2/software/miniconda/24.7.1/condabin:$PATH"
record="data/results/post_expansion/2026-10-21/reporting/executions/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$record"
exec > >(tee "$record/execution.log") 2>&1
trap 'code=$?; date -u > "$record/finished.txt"; printf "%s\n" "$code" > "$record/exit_code.txt"' EXIT
date -u > "$record/started.txt"
git rev-parse HEAD > "$record/git_head.txt"
git status --short > "$record/git_status.txt"
tar -czf "$record/source.tar.gz" config workflow scripts profiles envs
sha256sum "$record/source.tar.gz" > "$record/source_sha256.txt"
# Snakemake partitions the aggregation input list deterministically. Batches 1-6
# leave batch 7 and aggregation to the later unrestricted local evidence target.
python -m snakemake -s workflow/expanded_reporting.smk eukaryote_evidence \
  --batch "aggregate_ssu_read_lineages=$SLURM_ARRAY_TASK_ID/7" \
  --cores "$SLURM_CPUS_PER_TASK" --resources mem_mb=128000 --use-conda --conda-frontend mamba \
  --nolock --rerun-incomplete --rerun-triggers mtime --printshellcmds
