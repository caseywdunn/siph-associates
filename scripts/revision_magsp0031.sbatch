#!/usr/bin/env bash
#SBATCH --job-name=revision-MAGSP0031
#SBATCH --partition=ycga
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/revision-MAGSP0031-%j.log
set -euo pipefail
module reset
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/checkm2_env/bin:$PATH"
base="$PWD/data/results/deposit"
mag=Church2025__WAM-Z97910__bin.1
record="$base/executions/MAGSP0031-$SLURM_JOB_ID"
# Python multiprocessing uses Unix sockets with short pathname limits.
export TMPDIR="$PWD/tmp/c2-$SLURM_JOB_ID"
mkdir -p "$record" "$TMPDIR" "$base/cleaned_mags" "$base/edits"
git rev-parse HEAD > "$record/git_head.txt"
git diff > "$record/source.diff"
cp workflow/scripts/prepare_deposit_mag.py workflow/scripts/summarize_magsp0031_revision.py config/mag_contig_removals.tsv config/mag_adaptor_trims.tsv config/mag_sequence_repairs.tsv "$record/"
sha256sum "$record/"* > "$record/source.sha256"
trap 'status=$?; printf "%s\n" "$status" > "$record/exit_code.txt"' EXIT
python workflow/scripts/prepare_deposit_mag.py --mag "$mag" \
  --input "data/results/post_expansion/2026-10-20/curated/phase3_catalog/mags/candidates/$mag.fa" \
  --removals config/mag_contig_removals.tsv --trims config/mag_adaptor_trims.tsv \
  --repairs config/mag_sequence_repairs.tsv --output "$base/cleaned_mags/$mag.fa" \
  --report "$base/edits/$mag.json"
checkm2 --version > "$record/checkm2_version.txt" 2>&1
checkm2 predict --threads 8 --input "$base/cleaned_mags/$mag.fa" --extension fa \
  --database_path /gpfs/ycga/work/dunn/cwd7/databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd \
  --output-directory "$base/checkm2/$mag" --force

python workflow/scripts/summarize_magsp0031_revision.py --root "$base" --execution "$record"
