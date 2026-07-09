#!/usr/bin/env bash
# Exact read-pair count for one library: sum R1 reads across its fastq.gz files
# (R1 read count == read pairs for paired-end). Used for libraries not in SRA;
# SRA-deposited libraries use exact `spots` from runinfo instead.
# usage: count_reads.sh LIBRARY_ID "R1;R1;..." OUTDIR
set -euo pipefail
LIB="$1"; R1S="$2"; OUTDIR="$3"
mkdir -p "$OUTDIR"
source /etc/profile.d/z01_lmodinit.sh 2>/dev/null || true
module load pigz/2.7-GCCcore-12.2.0 2>/dev/null || true
DECOMP=$(command -v pigz >/dev/null && echo "pigz -dc -p ${SLURM_CPUS_PER_TASK:-4}" || echo "zcat")
total=0
IFS=';' read -ra FILES <<< "$R1S"
for f in "${FILES[@]}"; do
  [ -s "$f" ] || { echo "MISSING $f" >&2; exit 3; }
  n=$($DECOMP "$f" | wc -l)
  total=$((total + n))
done
pairs=$((total / 4))
printf '%s\t%s\n' "$LIB" "$pairs" > "$OUTDIR/${LIB//:/__}.count"
echo "[$(date -Iseconds)] $LIB read_pairs=$pairs (${#FILES[@]} files)"
