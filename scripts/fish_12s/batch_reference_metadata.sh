#!/usr/bin/env bash
#SBATCH --job-name=fish12s-taxonomy
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/fish_12s/taxonomy-%j.log
set -euo pipefail
module reset
refdir=data/results/fish_12s/reference/refseq_mito_2026-09-03
metadata=data/results/fish_12s/reference_metadata/2026-10-18
execution="$metadata/executions/$SLURM_JOB_ID"
mkdir -p "$execution"
cp scripts/fish_12s/{reference_metadata.py,batch_reference_metadata.sh} "$execution/"
date -u > "$execution/started.txt"
trap 'code=$?; date -u > "$execution/finished.txt"; printf "%s\n" "$code" > "$execution/exit_code.txt"' EXIT
if ! test -s "$refdir/mitochondrion.1.genomic.gbff.gz"; then
    curl --fail --location --retry 3 -D "$metadata/download.headers.txt" \
        https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.genomic.gbff.gz \
        -o "$refdir/mitochondrion.1.genomic.gbff.gz.part"
    gzip -t "$refdir/mitochondrion.1.genomic.gbff.gz.part"
    mv "$refdir/mitochondrion.1.genomic.gbff.gz.part" "$refdir/mitochondrion.1.genomic.gbff.gz"
fi
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python scripts/fish_12s/reference_metadata.py \
    --fasta "$refdir/mitochondrion.1.1.genomic.fna.gz" --genbank "$refdir/mitochondrion.1.genomic.gbff.gz" \
    --mitofish data/results/fish_12s/assessment/2026-10-17/reference_matches.tsv \
    --output "$metadata/references.tsv"
