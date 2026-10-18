#!/usr/bin/env bash
#SBATCH --job-name=fish12s-aligners
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --array=2,5
#SBATCH --output=logs/fish_12s/aligners.%A_%a.log
set -euo pipefail
module load Bowtie2/2.5.1-GCC-12.2.0
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python scripts/fish_12s/compare_aligners.py \
    --manifest "${1:-config/fish_12s_known_positive.tsv}" \
    --index "${SLURM_ARRAY_TASK_ID}" --threads 4 \
    --reference data/results/fish_12s/masked_reference/2026-10-17/fish_12s_masked.fasta \
    --bowtie-index data/results/fish_12s/positive_references/2026-10-17/fish_12s \
    --sortmerna data/results/fish_12s/software/sortmerna_env/bin/sortmerna \
    --output-root "${2:-data/results/fish_12s/positive_aligners/2026-10-17}"
