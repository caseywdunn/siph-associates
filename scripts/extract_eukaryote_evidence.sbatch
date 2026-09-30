#!/usr/bin/env bash
#SBATCH --job-name=sa.extract_eukaryote_evidence
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm/extract_eukaryote_evidence.%j.out
# Submit from the repository root: sbatch scripts/extract_eukaryote_evidence.sbatch
set -euo pipefail
python3 scripts/extract_eukaryote_evidence.py
