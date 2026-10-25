#!/usr/bin/env bash
#SBATCH --job-name=manuscript-revision-reporting
#SBATCH --partition=ycga
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=8:00:00
#SBATCH --output=logs/manuscript-revision-reporting-%j.log
set -euo pipefail
module reset
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:/vast/palmer/apps/avx2/software/miniconda/24.7.1/condabin:$PATH"
revision_root="data/results/manuscript_revision/2026-10-24"
export TMPDIR="$PWD/$revision_root/scratch"
export MPLCONFIGDIR="$PWD/$revision_root/matplotlib"
export XDG_CACHE_HOME="$PWD/$revision_root/cache"
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
record="$revision_root/executions/$SLURM_JOB_ID"
mkdir -p "$record"
exec > >(tee "$record/execution.log") 2>&1
trap 'code=$?; date -u > "$record/finished.txt"; printf "%s\n" "$code" > "$record/exit_code.txt"' EXIT
date -u > "$record/started.txt"
git rev-parse HEAD > "$record/git_head.txt"
git status --short > "$record/git_status.txt"
tar -czf "$record/source.tar.gz" config workflow scripts envs figures/*.py tests/test_contamination_policy.py
sha256sum "$record/source.tar.gz" > "$record/source_sha256.txt"
python tests/test_contamination_policy.py
python -m snakemake --snakefile workflow/manuscript_revision.smk revision_publication --dry-run --cores "$SLURM_CPUS_PER_TASK" --resources mem_mb=48000 --printshellcmds
for target in revision_cohort revision_statistics revision_publication; do
    python -m snakemake --snakefile workflow/manuscript_revision.smk "$target" \
        --cores "$SLURM_CPUS_PER_TASK" --resources mem_mb=48000 \
        --nolock --rerun-incomplete --rerun-triggers mtime --printshellcmds
    if [[ "$target" == revision_statistics ]]; then
        python scripts/export_domain_review_statistics.py
    fi
done
