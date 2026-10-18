#!/usr/bin/env bash
#SBATCH --job-name=fish12s-refs
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/fish_12s/refs.%j.log
set -euo pipefail
module load Bowtie2/2.5.1-GCC-12.2.0 BLAST+/2.15.0-gompi-2022b
refdir=${1:-data/results/fish_12s/positive_references/2026-10-17}
mkdir -p "$(dirname "$refdir")"
mkdir "$refdir"
cp scripts/fish_12s/prepare_positive_references.sbatch "$refdir/"
exec > >(tee "$refdir/preparation.log") 2>&1
set -x
date -u
bowtie2 --version
blastn -version
sha256sum data/results/fish_12s/reference/refseq_mito_2026-09-03/mitochondrion.1.1.genomic.fna.gz \
    data/results/fish_12s/masked_reference/2026-10-17/fish_12s_masked.fasta \
    data/results/fish_12s/assessment/2026-10-17/fish_12s.fasta > "$refdir/input_sha256.txt"
bowtie2-build --threads 4 \
    data/results/fish_12s/masked_reference/2026-10-17/fish_12s_masked.fasta "$refdir/fish_12s"
gzip -dc data/results/fish_12s/reference/refseq_mito_2026-09-03/mitochondrion.1.1.genomic.fna.gz \
    > "$refdir/competitive.fasta"
cat data/results/fish_12s/assessment/2026-10-17/fish_12s.fasta >> "$refdir/competitive.fasta"
makeblastdb -in "$refdir/competitive.fasta" -dbtype nucl -parse_seqids \
    -out "$refdir/competitive" -title 'RefSeq mitochondria 2026-09-03 and MitoFish 12S 2026-09'
sha256sum "$refdir/competitive.fasta" > "$refdir/output_sha256.txt"
date -u > "$refdir/COMPLETE"
