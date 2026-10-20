#!/usr/bin/env bash
#SBATCH --job-name=bargmannia-workflow
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/bargmannia/controller-%j.log
set -euo pipefail
module reset
module load BLAST+/2.15.0-gompi-2022b Bowtie2/2.5.1-GCC-12.2.0 SeqKit/2.8.1
export PATH=/gpfs/gibbs/project/dunn/cwd7/conda_envs/phyloflash_env/bin:$PATH
export PATH="$PATH:/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:/vast/palmer/apps/avx2/software/miniconda/24.7.1/condabin"
export XDG_CACHE_HOME="$PWD/.cache"
python=/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python
record="data/results/bargmannia/2026-10-19/executions/$SLURM_JOB_ID"
mkdir -p "$record"
exec > >(tee "$record/execution.log") 2>&1
trap 'code=$?; date -u > "$record/finished.txt"; printf "%s\n" "$code" > "$record/exit_code.txt"' EXIT
date -u > "$record/started.txt"
git rev-parse HEAD > "$record/git_head.txt"
git status --short > "$record/git_status.txt"
tar -czf "$record/source.tar.gz" Snakefile config workflow scripts profiles envs \
  manko2026_manifest.csv data/metadata/manko2026 data/sources/manko2026
sha256sum "$record/source.tar.gz" > "$record/source_sha256.txt"
"$python" scripts/fish_12s/check_cohort_controls.py \
  --root data/results/fish_12s/cohort_validation/2026-10-18 \
  --truth data/results/fish_12s/positive_aligner_controls/2026-10-17/synthetic/known_answer/checks.json \
  --reference-metadata data/results/fish_12s/reference_metadata/2026-10-18/references.tsv \
  --output "$record/FISH_CONTROLS_PASSED"
sha256sum -c data/results/fish_12s/positive_references/2026-10-17/output_sha256.txt
set -x
"$python" -m snakemake --snakefile workflow/bargmannia.smk --profile profiles/slurm \
  --jobs 4 --cores 10000 --local-cores 1 --nolock \
  --default-resources partition=ycga runtime=480 mem_mb=4000 \
  --resources cap_slots=1 --rerun-triggers mtime --rerun-incomplete
