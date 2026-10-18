#!/usr/bin/env bash
#SBATCH --job-name=coi-nr-audit
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/crustacean_coi/nr-audit-%j.log
set -euo pipefail
module reset
record=data/results/crustacean_coi/executions/${SLURM_JOB_ID}
output=data/results/crustacean_coi/validation/2026-10-17/nr_fast_audit
mkdir -p "$record" "$output"
exec > >(tee "$record/execution.log") 2>&1
trap 'code=$?; date -u > "$record/finished.txt"; printf "%s\n" "$code" > "$record/exit_code.txt"' EXIT
cp scripts/crustacean_coi/batch_nr_audit.sh "$record/"
date -u > "$record/started.txt"
git rev-parse HEAD > "$record/git_head.txt"
stat /gpfs/ycga/work/dunn/cwd7/databases/diamond/nr.dmnd > "$record/database_stat.txt"
sha256sum data/results/crustacean_coi/validation/2026-10-17/audit/queries.fna > "$record/input_sha256.txt"
test ! -e "$output/hits.tsv"
set -x
/vast/palmer/apps/avx2/software/SeqKit/2.8.1/bin/seqkit grep --threads 1 -r -p 'NA35__c000000(08|16)_[12]$|^REF__' \
    data/results/crustacean_coi/validation/2026-10-17/audit/queries.fna > "$output/queries.fna"
/vast/palmer/apps/avx2/software/DIAMOND/2.1.24-GCC-13.3.0/bin/diamond version
/vast/palmer/apps/avx2/software/DIAMOND/2.1.24-GCC-13.3.0/bin/diamond dbinfo --db /gpfs/ycga/work/dunn/cwd7/databases/diamond/nr.dmnd
/vast/palmer/apps/avx2/software/DIAMOND/2.1.24-GCC-13.3.0/bin/diamond blastx \
    --query "$output/queries.fna" --db /gpfs/ycga/work/dunn/cwd7/databases/diamond/nr.dmnd \
    --out "$output/hits.tsv" --query-gencode 5 --fast --min-orf 20 --evalue 1e-3 \
    --query-cover 60 --max-target-seqs 25 --block-size 4 --index-chunks 1 --threads "$SLURM_CPUS_PER_TASK" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle
sha256sum "$output/queries.fna" "$output/hits.tsv" > "$record/output_sha256.txt"
