#!/usr/bin/env bash
#SBATCH --job-name=sa.extract_eukaryote_read_evidence
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/extract_eukaryote_read_evidence.%j.out
# Submit from the repository root: sbatch scripts/extract_eukaryote_read_evidence.sbatch
set -euo pipefail
python3 scripts/extract_eukaryote_read_evidence.py
