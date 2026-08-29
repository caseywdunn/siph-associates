#!/usr/bin/env bash
set -euo pipefail

PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SNAKEMAKE=/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/snakemake
export XDG_CACHE_HOME="$PROJECT/.cache"
mkdir -p "$XDG_CACHE_HOME"
cd "$PROJECT"

python3 scripts/build_workflow_samples.py
python3 scripts/build_fixture_samples.py
python3 scripts/validate_manifest.py --freeze data/metadata/manifest.freeze.sha256
"$SNAKEMAKE" --lint --configfile config/config.yaml
"$SNAKEMAKE" phase1_ready -n --cores 1 --configfile config/config.yaml
"$SNAKEMAKE" screen_cohort -n --quiet --cores 1 --configfile config/config.yaml
"$SNAKEMAKE" phase1_smoke -n --quiet --cores 1 --configfile tests/config.fixture.yaml
python3 scripts/validate_phase1.py "$@"
