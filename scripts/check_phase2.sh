#!/usr/bin/env bash
set -euo pipefail

PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SNAKEMAKE=/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/snakemake
export XDG_CACHE_HOME="$PROJECT/.cache"
mkdir -p "$XDG_CACHE_HOME"
cd "$PROJECT"

python3 scripts/build_workflow_samples.py
python3 scripts/validate_manifest.py --freeze data/metadata/manifest.freeze.sha256
"$SNAKEMAKE" --lint --configfile config/config.yaml
"$SNAKEMAKE" screen_pilot -n --quiet --profile profiles/slurm --configfile config/config.yaml
"$SNAKEMAKE" screen_cohort -n --quiet --profile profiles/slurm --configfile config/config.yaml
"$SNAKEMAKE" \
  "$PROJECT/tests/cap_scratch/trimmed/Church2025__FM-16644_R1.fastq.gz" \
  "$PROJECT/tests/cap_scratch/trimmed/Church2025__FM-16644_R2.fastq.gz" \
  --cores 6 --use-envmodules --rerun-incomplete -R trim_reads \
  --configfile tests/config.cap.fixture.yaml
python3 tests/check_capped_fixture.py
python3 tests/test_atomic_move.py
python3 tests/test_phase2_fixture.py
python3 scripts/validate_phase2.py "$@"
