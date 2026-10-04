#!/usr/bin/env bash
# Bounded publication-review rerun: accepted sequence outputs remain inputs.
set -euo pipefail
PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$PROJECT"
RULES='normalize_phase6_metadata test_contamination_flowcell test_nanomia_contamination_flowcell summarize_incidence fit_physalia_models_lme4 analyze_presence_threshold compare_sensitivity_analyses validate_phase6 phase6_analysis grade_eukaryotes validate_eukaryote_gate eukaryote_gate publication_cohort_figures publication_mycoplasmatales_figures'
bash scripts/submit_workflow.sh \
    'phase6_analysis eukaryote_gate publication_cohort_figures publication_mycoplasmatales_figures' \
    config/config.yaml "$RULES"
