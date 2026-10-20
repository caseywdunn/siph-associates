#!/usr/bin/env bash
#SBATCH --job-name=mag-domain-audit
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/mag_domain_audit/audit-%j.log
set -euo pipefail
python=/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python
run=$("$python" -c 'import json; print(json.load(open("config/mag_domain_audit.json"))["run"])')
execution="$run/executions/$SLURM_JOB_ID"
mkdir -p "$execution/source"
date -u > "$execution/started.txt"
git rev-parse HEAD > "$execution/git_head.txt"
git status --short > "$execution/git_status.txt"
tar -czf "$execution/source.tar.gz" scripts/mag_domain_audit workflow/mag_domain_audit.smk config/mag_domain_audit.json
sha256sum "$execution/source.tar.gz" > "$execution/source_sha256.txt"
tar -xzf "$execution/source.tar.gz" -C "$execution/source"
exec > >(tee "$execution/execution.log") 2>&1
trap 'code=$?; date -u > "$execution/finished.txt"; printf "%s\n" "$code" > "$execution/exit_code.txt"' EXIT
cd "$execution/source"
export XDG_CACHE_HOME="$execution/cache"
"$python" - <<'PY'
import json, pathlib, subprocess
config = json.loads(pathlib.Path("config/mag_domain_audit.json").read_text())
for executable, argument in [(config['diamond'], 'version'), (config['gtdbtk_prefix']+'/bin/gtdbtk', '--version')]:
    subprocess.run([executable, argument], check=True)
prefixes = {pathlib.Path(config['python']).parent.parent, pathlib.Path(config['diamond']).parent.parent, pathlib.Path(config['gtdbtk_prefix'])}
records = {str(prefix): [json.loads(p.read_text()) for p in sorted((prefix/'conda-meta').glob('*.json'))] for prefix in prefixes}
pathlib.Path('../conda_packages.json').write_text(json.dumps(records, indent=2)+'\n')
PY
set -x
"$python" -m snakemake all --snakefile workflow/mag_domain_audit.smk \
    --configfile config/mag_domain_audit.json --cores "$SLURM_CPUS_PER_TASK" \
    --resources mem_mb=60000 --printshellcmds --rerun-incomplete --latency-wait 60
find "$run/review" "$run/controls" -type f -print0 | sort -z | xargs -0 sha256sum > "$execution/output_sha256.txt"
