#!/usr/bin/env bash
#SBATCH --job-name=sa.phase3_qc_evidence
#SBATCH --partition=day
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=02:00:00
#SBATCH --output=logs/slurm/phase3_qc_evidence.%j.out
# Regenerate the evidence for docs/phase3_qc_decisions.md from accepted outputs.
# Submit from the repository root: sbatch scripts/phase3_qc_evidence.sbatch
set -euo pipefail
BIN=/gpfs/gibbs/project/dunn/cwd7/conda_envs/assembly_env/bin
N_SEPTATA=/gpfs/ycga/work/dunn/cwd7/databases/refgenomes/nanomia_septata/GCA_048301705.1_ASM4830170v1_genomic.fna
P_PHYSALIS=/gpfs/ycga/work/dunn/cwd7/databases/refgenomes/physalia_physalis/GCA_041430235.2_Physalia_TX2017-38_primary_02_genomic.fna
OUT=data/results/phase3_cohort/qc_evidence

python3 scripts/phase3_qc_evidence.py prepare
# Host-genome test: do viral contigs from reference-free Nanomia libraries align to the congener reference?
for query in nanomia_free other_free septata_depleted; do
  "$BIN/minimap2" -x asm20 -t 8 --secondary=no -c "$N_SEPTATA" "$OUT/$query.fna" > "$OUT/$query.vs_septata.paf"
done
"$BIN/minimap2" -x asm20 -t 8 --secondary=no -c "$P_PHYSALIS" "$OUT/physalia_5kb.fa" > "$OUT/physalia_5kb.vs_physalis.paf"
# Pairwise ANI among medium-quality-or-better MAGs.
"$BIN/skani" triangle -t 8 -E "$OUT"/medium_bins/*.fa > "$OUT/medium_bins_skani.tsv"
"$BIN/minimap2" --version
"$BIN/skani" --version
python3 scripts/phase3_qc_evidence.py summarize
