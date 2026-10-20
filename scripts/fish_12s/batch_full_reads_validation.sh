#!/usr/bin/env bash
#SBATCH --job-name=fish12s-full-controls
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=logs/fish_12s/full-controls-%j.log
set -euo pipefail
python=/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python
"$python" tests/integration/fish_full_reads_fixture.py prepare
bash scripts/fish_12s/batch_full_reads.sh data/results/fish_12s/full_reads_validation/2026-10-18/capped.json
bash scripts/fish_12s/batch_full_reads.sh data/results/fish_12s/full_reads_validation/2026-10-18/uncapped.json
"$python" tests/integration/fish_full_reads_fixture.py check
