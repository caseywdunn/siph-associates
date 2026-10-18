#!/usr/bin/env bash
#SBATCH --job-name=fish12s-positive
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-8%4
#SBATCH --output=logs/fish_12s/positive.%A_%a.log
set -euo pipefail
module load SeqKit/2.8.1
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/phyloflash_env/bin:$PATH"
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python scripts/fish_12s/screen_positive.py \
    --manifest "${1:-config/fish_12s_known_positive.tsv}" \
    --index "${SLURM_ARRAY_TASK_ID}" --threads 4 \
    --reference data/results/fish_12s/masked_reference/2026-10-17/fish_12s_masked.fasta \
    --output-root "${2:-data/results/fish_12s/known_positive/2026-10-17}"
