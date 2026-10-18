#!/usr/bin/env bash
#SBATCH --job-name=fish12s-verify
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --array=0-8%4
#SBATCH --output=logs/fish_12s/verify.%A_%a.log
set -euo pipefail
module load BLAST+/2.15.0-gompi-2022b
manifest=${1:-config/fish_12s_known_positive.tsv}
positive_root=${2:-data/results/fish_12s/known_positive/2026-10-17}
read -r sample lane < <(/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python - "$manifest" "$SLURM_ARRAY_TASK_ID" <<'PY'
import csv
import sys
with open(sys.argv[1]) as handle:
    row = list(csv.DictReader(handle, delimiter='\t'))[int(sys.argv[2])]
print(row['analysis_sample_id'], row['read_pair_id'])
PY
)
input_dir="$positive_root/$sample/$lane"
out="$input_dir/alignment"
mkdir "$out"
cp scripts/fish_12s/align_positive_candidates.sbatch "$out/"
exec > >(tee "$out/run.log") 2>&1
set -x
date -u
blastn -version
sha256sum "$input_dir/candidates.fasta" "$input_dir/run.json" \
    data/results/fish_12s/positive_references/2026-10-17/input_sha256.txt \
    data/results/fish_12s/positive_references/2026-10-17/output_sha256.txt > "$out/input_sha256.txt"
if test -s "$input_dir/candidates.fasta"; then
    /usr/bin/time -v blastn -task blastn -word_size 11 -dust yes -soft_masking true \
        -query "$input_dir/candidates.fasta" \
        -db data/results/fish_12s/positive_references/2026-10-17/competitive \
        -evalue 1e-8 -max_target_seqs 200 -max_hsps 1 -num_threads 4 \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle' \
        -out "$out/hits.tsv"
else
    touch "$out/hits.tsv"
fi
sha256sum "$out/hits.tsv" > "$out/output_sha256.txt"
date -u > "$out/COMPLETE"
