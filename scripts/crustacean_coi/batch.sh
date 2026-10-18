#!/usr/bin/env bash
#SBATCH --job-name=crustacean-coi
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --output=logs/crustacean_coi/%j.log
set -euo pipefail
module reset
# Load BLAST dependencies. DIAMOND and SeqKit are invoked by absolute path.
module load BLAST+/2.15.0-gompi-2022b
export PATH=/gpfs/gibbs/project/dunn/cwd7/conda_envs/phyloflash_env/bin:$PATH
target=${1:-pilot}
run_record=data/results/crustacean_coi/executions/${SLURM_JOB_ID}
mkdir -p "$run_record"
exec > >(tee "$run_record/execution.log") 2>&1
trap 'code=$?; date -u > "$run_record/finished.txt"; printf "%s\n" "$code" > "$run_record/exit_code.txt"' EXIT
date -u > "$run_record/started.txt"
git rev-parse HEAD > "$run_record/git_head.txt"
git status --short > "$run_record/git_status.txt"
tar -czf "$run_record/source.tar.gz" scripts/crustacean_coi workflow/crustacean_coi.smk workflow/rules/crustacean_coi_search.smk config/crustacean_coi*
sha256sum "$run_record/source.tar.gz" > "$run_record/source.sha256"
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python - "$run_record" <<'PY'
import csv
import json
import sys
from pathlib import Path
import Bio
import snakemake
records = []
with open('config/crustacean_coi_pilot.tsv') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        for mate in ('r1', 'r2'):
            path = Path(row[f'{mate}_path'])
            stat = path.stat()
            if stat.st_size != int(row[f'{mate}_bytes']):
                raise ValueError(f'Raw input size changed: {path}')
            records.append({'sample': row['analysis_sample_id'], 'mate': mate,
                            'path': str(path), 'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns})
Path(sys.argv[1], 'inputs.json').write_text(json.dumps(records, indent=2) + '\n')
Path(sys.argv[1], 'python_versions.json').write_text(json.dumps({
    'python': sys.version, 'biopython': Bio.__version__, 'snakemake': snakemake.__version__
}, indent=2) + '\n')
PY
set -x
/vast/palmer/apps/avx2/software/DIAMOND/2.1.24-GCC-13.3.0/bin/diamond version
blastn -version
/vast/palmer/apps/avx2/software/SeqKit/2.8.1/bin/seqkit version
bbduk.sh --version
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python -m snakemake \
    --snakefile workflow/crustacean_coi.smk "$target" \
    --cores "${SLURM_CPUS_PER_TASK}" --resources mem_mb=18000 \
    --printshellcmds --rerun-incomplete --rerun-triggers mtime --latency-wait 60
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python - "$run_record" <<'PY'
import hashlib
import json
import shutil
import sys
from pathlib import Path
configuration = json.loads(Path('config/crustacean_coi.json').read_text())
destination = Path(sys.argv[1], 'reports')
for key in ('panel', 'run'):
    root = Path(configuration[key])
    for path in root.rglob('*'):
        if path.is_file() and (path.suffix == '.json' or path.name in (
            'coverage.tsv', 'classified_pairs.tsv', 'candidate_pairs.tsv', 'PASSED'
        ) or path.parent.name == 'benchmarks'):
            target = destination / key / path.relative_to(root)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)
with Path(sys.argv[1], 'output_sha256.txt').open('w') as handle:
    for path in sorted(destination.rglob('*')):
        if path.is_file():
            with path.open('rb') as stream:
                handle.write(f'{hashlib.file_digest(stream, "sha256").hexdigest()}  {path}\n')
PY
