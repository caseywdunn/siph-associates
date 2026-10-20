#!/usr/bin/env bash
#SBATCH --job-name=catalog-domain-review
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/catalog_domain_review/review-%j.log
set -euo pipefail
export XDG_CACHE_HOME="$PWD/.cache"
export PATH=/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:$PATH
record="data/results/catalog_domain_review/2026-10-19/executions/$SLURM_JOB_ID"
mkdir -p "$record"
exec > >(tee "$record/execution.log") 2>&1
trap 'code=$?; date -u > "$record/finished.txt"; printf "%s\n" "$code" > "$record/exit_code.txt"' EXIT
date -u > "$record/started.txt"
git rev-parse HEAD > "$record/git_head.txt"
git status --short > "$record/git_status.txt"
tar -czf "$record/source.tar.gz" Snakefile config workflow scripts figures/*.py
sha256sum "$record/source.tar.gz" data/results/catalog_domain_review/2026-10-19/pre_correction.tar.gz > "$record/source_sha256.txt"
set -x
snakemake phase6_analysis publication_cohort_figures publication_collection_context_figures \
  --cores 12 --resources mem_mb=46000 --rerun-triggers mtime --nolock \
  --printshellcmds --rerun-incomplete \
  --allowed-rules grade_presences test_contamination_flowcell test_nanomia_contamination_flowcell \
    summarize_incidence fit_physalia_models_lme4 compare_sensitivity_analyses \
    analyze_presence_threshold validate_phase6 phase6_analysis apply_catalog_domain_review \
    publication_cohort_figures collection_context_tables publication_collection_context_figures
python figures/plot_associate_overview.py
python scripts/summarize_dalila_revision.py
python scripts/summarize_screening_results.py
python scripts/summarize_mdb_revision.py
python scripts/export_domain_review_statistics.py
sha256sum figures/cohort/summary.json figures/catalog_review/summary.json \
  figures/screening/summary.json figures/collection_context/summary.json \
  data/results/stages/phase6_analysis_v2.done > "$record/output_sha256.txt"
