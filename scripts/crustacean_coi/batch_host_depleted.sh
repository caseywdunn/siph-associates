#!/usr/bin/env bash
#SBATCH --job-name=coi-host-depleted
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/crustacean_coi/host-depleted-%j.log
set -euo pipefail
module reset
module load BLAST+/2.15.0-gompi-2022b
export PATH=/gpfs/gibbs/project/dunn/cwd7/conda_envs/phyloflash_env/bin:$PATH
run_record=data/results/crustacean_coi/executions/${SLURM_JOB_ID}
mkdir -p "$run_record"
exec > >(tee "$run_record/execution.log") 2>&1
trap 'code=$?; date -u > "$run_record/finished.txt"; printf "%s\n" "$code" > "$run_record/exit_code.txt"' EXIT
date -u > "$run_record/started.txt"
pwd > "$run_record/working_directory.txt"
git rev-parse HEAD > "$run_record/git_head.txt"
git status --short > "$run_record/git_status.txt"
tar -czf "$run_record/source.tar.gz" scripts/crustacean_coi workflow/crustacean_coi*.smk workflow/rules/crustacean_coi*.smk config/crustacean_coi*.json config/crustacean_coi_pilot.tsv config/prey_marker_host_depleted.tsv
sha256sum "$run_record/source.tar.gz" > "$run_record/source.sha256"
sha256sum data/results/crustacean_coi/panel/2026-10-17/{target.masked.fna,target.dmnd,competitive.nsq,competitive.nhr,competitive.nin} \
    data/results/crustacean_coi/verification/2026-10-18/{gate.json,curated/full.dmnd,curated/reference_metadata.tsv,curated/protein_metadata.tsv} \
    data/results/crustacean_coi/host_depleted/2026-10-17/Ahuja2026__*/{candidates.fna,candidate_pairs.tsv,competitive.tsv,recruitment.json} > "$run_record/reference_and_reused_input_sha256.txt"
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python - "$run_record" <<'PY'
import csv,json,sys
from pathlib import Path
with open('config/prey_marker_host_depleted.tsv') as handle:
    rows=list(csv.DictReader(handle,delimiter='\t'))
records={r[f'r{m}_path']: dict(size=Path(r[f'r{m}_path']).stat().st_size,mtime_ns=Path(r[f'r{m}_path']).stat().st_mtime_ns) for r in rows for m in [1,2]}
Path(sys.argv[1],'input_stats.json').write_text(json.dumps(records,indent=2)+'\n')
PY
set -x
/vast/palmer/apps/avx2/software/DIAMOND/2.1.24-GCC-13.3.0/bin/diamond version
blastn -version
bbduk.sh --version
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python -m snakemake host_depleted \
    --snakefile workflow/crustacean_coi_host_depleted.smk --cores "$SLURM_CPUS_PER_TASK" \
    --resources mem_mb=90000 --printshellcmds --rerun-incomplete --rerun-triggers mtime --latency-wait 60
rg --files -0 data/results/crustacean_coi/cohort_host_depleted/2026-10-18/verified | sort -z | xargs -0 sha256sum > "$run_record/verified_output_sha256.txt"
