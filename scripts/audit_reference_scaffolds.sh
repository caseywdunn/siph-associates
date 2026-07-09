#!/usr/bin/env bash
#SBATCH --job-name=ref_audit
#SBATCH --partition=day
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/ref_audit_%j.out
# S3a reference-scaffold audit: extract a reference's unplaced scaffolds and run
# geNomad end-to-end to confirm no co-assembled bacterial/symbiont sequence sits
# in the host reference (which would wrongly delete real symbiont reads under
# whole-genome host subtraction), and to report co-assembled non-host / host EVEs.
#   usage: sbatch scripts/audit_reference_scaffolds.sh GENOME.fna.gz UNPLACED_NAMES.txt OUTDIR
set -euo pipefail
GENOME_GZ="$1"; UNPLACED="$2"; OUTDIR="$3"
ENV=/gpfs/gibbs/project/dunn/cwd7/conda_envs/genomad_env
AENV=/gpfs/gibbs/project/dunn/cwd7/conda_envs/assembly_env
GENOMAD_DB=/gpfs/ycga/work/dunn/cwd7/databases/genomad_db/genomad_db
mkdir -p "$OUTDIR"
FNA="${GENOME_GZ%.gz}"; UNP="$OUTDIR/unplaced.fasta"
export PATH="$AENV/bin:$PATH"
[[ -s "$FNA" ]] || gunzip -kf "$GENOME_GZ"
samtools faidx "$FNA"
echo "[$(date -Iseconds)] extracting $(wc -l < "$UNPLACED") unplaced scaffolds"
samtools faidx "$FNA" $(cat "$UNPLACED") > "$UNP"
samtools faidx "$UNP"
echo "[$(date -Iseconds)] unplaced total: $(awk '{s+=$2}END{printf "%.1f Mb\n",s/1e6}' "$UNP.fai")"
source /etc/profile.d/z01_lmodinit.sh 2>/dev/null || true
module load miniconda/24.9.2
set +u; conda activate "$ENV"; set -u
echo "[$(date -Iseconds)] geNomad $(genomad --version 2>&1 | head -1)"
genomad end-to-end --threads "${SLURM_CPUS_PER_TASK:-16}" --cleanup "$UNP" "$OUTDIR/genomad" "$GENOMAD_DB"
echo "[$(date -Iseconds)] done -> $OUTDIR/genomad"
