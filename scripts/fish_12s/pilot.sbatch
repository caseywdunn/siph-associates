#!/usr/bin/env bash
#SBATCH --job-name=fish12s-pilot
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --array=0-3%2
#SBATCH --output=logs/fish_12s/pilot.%A_%a.log
set -euo pipefail
module load SeqKit/2.8.1
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/phyloflash_env/bin:$PATH"
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python scripts/fish_12s/pilot.py \
    --manifest "${1:-config/fish_12s_pilot.tsv}" \
    --index "${SLURM_ARRAY_TASK_ID}" --pairs 1000000 \
    --reference data/results/fish_12s/assessment/2026-10-17/fish_12s.fasta \
    --output-root data/results/fish_12s/pilot/2026-10-17 \
    --bbduk /gpfs/gibbs/project/dunn/cwd7/conda_envs/phyloflash_env/bin/bbduk.sh
