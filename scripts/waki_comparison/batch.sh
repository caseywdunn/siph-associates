#!/usr/bin/env bash
#SBATCH --job-name=sa.waki_comparison
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/slurm/waki_comparison.%j.out
set -euo pipefail
module reset
cd "${SIPH_ASSOCIATES_PROJECT:-/gpfs/ycga/work/dunn/cwd7/siph_associates}"
WAKI_RUN_DIR="${WAKI_RUN_DIR:-data/results/waki_comparison/reproduction-${SLURM_JOB_ID}}"
WAKI_PYTHON=/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python
"$WAKI_PYTHON" scripts/waki_comparison/compare_waki_assemblies.py \
  --include-context --output "$WAKI_RUN_DIR"
"$WAKI_PYTHON" scripts/waki_comparison/summarize_waki_comparison.py \
  --analysis "$WAKI_RUN_DIR" --output "$WAKI_RUN_DIR/report"
