#!/usr/bin/env bash
#SBATCH --job-name=cohort-remapping-v3
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/cohort-remapping-%j.log
set -euo pipefail
module reset
export XDG_CACHE_HOME="$PWD/.cache"
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:/vast/palmer/apps/avx2/software/miniconda/24.7.1/condabin:$PATH"
record="data/results/post_expansion/2026-10-20/remapping/executions/$SLURM_JOB_ID"
mkdir -p "$record"
exec > >(tee "$record/execution.log") 2>&1
trap 'code=$?; date -u > "$record/finished.txt"; printf "%s\n" "$code" > "$record/exit_code.txt"' EXIT
date -u > "$record/started.txt"
git rev-parse HEAD > "$record/git_head.txt"
git status --short > "$record/git_status.txt"
tar -czf "$record/source.tar.gz" config workflow scripts profiles envs
sha256sum "$record/source.tar.gz" > "$record/source_sha256.txt"
python - "$record/mapping_software.json" <<'PY'
import json
import subprocess
import sys
from pathlib import Path

prefix = Path('/gpfs/gibbs/project/dunn/cwd7/conda_envs/assembly_env/bin')
versions = {}
for name, flags, expected in [('bwa', [], 1), ('samtools', ['--version'], 0), ('coverm', ['--version'], 0)]:
    result = subprocess.run([str(prefix / name), *flags], capture_output=True, text=True, check=False)
    if result.returncode != expected:
        raise RuntimeError(f'{name} version check failed: {result.stderr}')
    versions[name] = (result.stdout + result.stderr).strip()
Path(sys.argv[1]).write_text(json.dumps(versions, indent=2) + '\n')
PY
set -x
python -m snakemake --snakefile workflow/cohort_remapping.smk "${1:-remapping}" \
  --profile profiles/slurm --jobs 16 --cores 256 --resources cap_slots=1 \
  --nolock --rerun-triggers mtime
