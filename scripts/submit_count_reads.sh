#!/usr/bin/env bash
#SBATCH --job-name=count_reads
#SBATCH --partition=day
#SBATCH --cpus-per-task=6
#SBATCH --mem=4G
#SBATCH --time=08:00:00
#SBATCH --array=1-21%8
#SBATCH --output=logs/count_reads_%A_%a.out
# Exact read-pair counts for the libraries not in SRA (data/sources/reads_to_count.tsv).
#   usage: sbatch scripts/submit_count_reads.sh
set -euo pipefail
PROJECT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
LIST="$PROJECT/data/sources/reads_to_count.tsv"
OUTDIR="$PROJECT/data/results/read_counts"
ROW=$(tail -n +2 "$LIST" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[[ -n "${ROW:-}" ]] || { echo "no library at task $SLURM_ARRAY_TASK_ID" >&2; exit 3; }
LIB=$(cut -f1 <<<"$ROW"); R1S=$(cut -f3 <<<"$ROW")
echo "[$(date -Iseconds)] task=$SLURM_ARRAY_TASK_ID $LIB"
bash "$PROJECT/scripts/count_reads.sh" "$LIB" "$R1S" "$OUTDIR"
