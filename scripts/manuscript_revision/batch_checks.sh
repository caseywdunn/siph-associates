#!/usr/bin/env bash
#SBATCH --job-name=manuscript-revision-checks
#SBATCH --partition=ycga
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=logs/manuscript-revision-checks-%j.log
set -euo pipefail
module reset
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:/vast/palmer/apps/avx2/software/miniconda/24.7.1/condabin:$PATH"
record="data/results/manuscript_revision/2026-10-24/checks/$SLURM_JOB_ID"
mkdir -p "$record/tmp"
export TMPDIR="$PWD/$record/tmp"
export MPLCONFIGDIR="$PWD/$record/matplotlib"
export XDG_CACHE_HOME="$PWD/$record/cache"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
exec > >(tee "$record/execution.log") 2>&1
trap 'code=$?; date -u > "$record/finished.txt"; printf "%s\n" "$code" > "$record/exit_code.txt"' EXIT
date -u > "$record/started.txt"
git rev-parse HEAD > "$record/git_head.txt"
git status --short > "$record/git_status.txt"
tar -czf "$record/source.tar.gz" config workflow scripts tests figures/*.py
sha256sum "$record/source.tar.gz" > "$record/source_sha256.txt"
python -m pytest -q --basetemp="$TMPDIR/pytest" -o cache_dir="$PWD/$record/pytest_cache" \
  --junitxml="$record/tests.xml" \
  tests/test_phase6_statistics.py tests/test_eukaryote_reporting.py \
  tests/test_viral_reporting_classification.py tests/test_coauthor_metadata.py \
  tests/test_collection_context.py tests/test_expanded_reporting.py \
  tests/test_reporting_manifest.py tests/test_contamination_policy.py \
  tests/test_fish_12s_reporting.py tests/test_human_cutibacterium.py
python - "$record" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

manuscript_scripts = Path('/gpfs/ycga/work/dunn/cwd7/manuscript_siph_associates/scripts')
sys.path.insert(0, str(manuscript_scripts))
import sync_assets

result = {'source': str(manuscript_scripts/'sync_assets.py'),
          'source_sha256': hashlib.sha256((manuscript_scripts/'sync_assets.py').read_bytes()).hexdigest(),
          'checks': {}}
for name in ('recorded_figure_hashes', 'recorded_screening_hashes'):
    try:
        hashes = getattr(sync_assets, name)(Path.cwd())
        result['checks'][name] = {'status': 'PASS', 'verified_files': len(hashes)}
    except Exception as error:
        result['checks'][name] = {'status': 'FAIL', 'error': str(error)}
result['status'] = 'PASS' if all(v['status']=='PASS' for v in result['checks'].values()) else 'FAIL'
Path(sys.argv[1], 'provenance_checks.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
if result['status'] != 'PASS':
    raise SystemExit(1)
PY
