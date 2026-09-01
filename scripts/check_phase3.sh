#!/usr/bin/env bash
set -euo pipefail

PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
SNAKEMAKE=/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/snakemake
export XDG_CACHE_HOME="$PROJECT/.cache"
export CONDA_PKGS_DIRS="$PROJECT/.cache/conda/pkgs"
export PATH="/vast/palmer/apps/avx2/software/miniconda/24.7.1/condabin:$PATH"
mkdir -p "$XDG_CACHE_HOME" "$CONDA_PKGS_DIRS"
cd "$PROJECT"

python3 scripts/validate_manifest.py --freeze data/metadata/manifest.freeze.sha256
python3 scripts/validate_phase2.py --scope cohort
python3 scripts/validate_phase3.py "$@"
"$SNAKEMAKE" --lint --configfile config/config.yaml
"$SNAKEMAKE" phase3_benchmark -n --quiet --profile profiles/slurm --configfile config/config.yaml
