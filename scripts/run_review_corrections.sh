#!/usr/bin/env bash
# Run the bounded workflow inside an existing compute allocation (no submission).
set -euo pipefail
: "${SLURM_JOB_ID:?Run inside a SLURM compute allocation}"
PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$PROJECT"
export XDG_CACHE_HOME="$PROJECT/.cache"
export CONDA_PKGS_DIRS="$PROJECT/.cache/conda/pkgs"
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:/vast/palmer/apps/avx2/software/miniconda/24.7.1/condabin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p logs/review "$XDG_CACHE_HOME" "$CONDA_PKGS_DIRS"
# Sequential execution fits the checked allocation (one CPU and 5 GB).
exec snakemake --configfile config/config.yaml --cores 1 --local-cores 1 \
    --resources mem_mb=4096 --use-conda --conda-frontend mamba \
    --rerun-incomplete --rerun-triggers mtime --latency-wait 60 --printshellcmds \
    "$@" \
    --allowed-rules normalize_phase6_metadata test_contamination_flowcell \
    test_nanomia_contamination_flowcell summarize_incidence fit_physalia_models_lme4 \
    analyze_presence_threshold compare_sensitivity_analyses validate_phase6 phase6_analysis \
    grade_eukaryotes validate_eukaryote_gate eukaryote_gate \
    publication_cohort_figures publication_mycoplasmatales_figures \
    -- phase6_analysis eukaryote_gate publication_cohort_figures publication_mycoplasmatales_figures
