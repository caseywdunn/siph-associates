#!/usr/bin/env bash
#SBATCH --job-name=fish12s-host-depleted
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=72G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/fish_12s/host-depleted-%j.log
set -euo pipefail
module reset
module load BLAST+/2.15.0-gompi-2022b Bowtie2/2.5.1-GCC-12.2.0 SeqKit/2.8.1
export PATH=/gpfs/gibbs/project/dunn/cwd7/conda_envs/phyloflash_env/bin:$PATH
run_record=data/results/fish_12s/executions/${SLURM_JOB_ID}
fish_config=${1:-config/fish_12s_host_depleted_cohort.json}
fish_target=${2:-host_depleted}
fish_root=$(/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python -c 'import json,sys; print(json.load(open(sys.argv[1]))["run"])' "$fish_config")
fish_manifest=$(/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python -c 'import json,sys; print(json.load(open(sys.argv[1]))["manifest"])' "$fish_config")
mkdir -p "$run_record"
exec > >(tee "$run_record/execution.log") 2>&1
trap 'code=$?; date -u > "$run_record/finished.txt"; printf "%s\n" "$code" > "$run_record/exit_code.txt"' EXIT
date -u > "$run_record/started.txt"
pwd > "$run_record/working_directory.txt"
git rev-parse HEAD > "$run_record/git_head.txt"
git status --short > "$run_record/git_status.txt"
tar -czf "$run_record/source.tar.gz" scripts/fish_12s workflow/fish_12s_host_depleted.smk workflow/rules/fish_12s_search.smk config/fish_12s_host_depleted_cohort.json config/fish_12s_host_depleted.tsv config/prey_marker_host_depleted.tsv config/fish_12s_cohort_validation.* "$fish_config" "$fish_manifest"
sha256sum "$run_record/source.tar.gz" > "$run_record/source.sha256"
sha256sum data/results/fish_12s/positive_references/2026-10-17/{input_sha256.txt,output_sha256.txt} \
    data/results/fish_12s/masked_reference/2026-10-17/fish_12s_masked.fasta > "$run_record/reference_sha256.txt"
sha256sum -c data/results/fish_12s/positive_references/2026-10-17/output_sha256.txt
# Production screening requires the fish workflow's own completed controls.
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python - "$fish_config" <<'PY'
import json
import csv
import sys
import subprocess
from pathlib import Path
config=json.loads(Path(sys.argv[1]).read_text())
if config.get('synthetic_validation'):
    with Path(config['manifest']).open() as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    if len(rows) != 1 or rows[0]['analysis_sample_id'] != 'synthetic' or int(rows[0]['pairs']) != 9:
        raise SystemExit('Unexpected synthetic validation input')
    sys.exit(0)
subprocess.run([sys.executable, 'scripts/fish_12s/check_cohort_controls.py',
    '--root', 'data/results/fish_12s/cohort_validation/2026-10-18',
    '--truth', 'data/results/fish_12s/positive_aligner_controls/2026-10-17/synthetic/known_answer/checks.json',
    '--reference-metadata', 'data/results/fish_12s/reference_metadata/2026-10-18/references.tsv',
    '--output', 'data/results/fish_12s/control_review/2026-10-18/CONTROLS_PASSED'], check=True)
PY
set -x
blastn -version
bowtie2 --version
seqkit version
bbduk.sh --version
data/results/fish_12s/software/sortmerna_env/bin/sortmerna --version
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python -m snakemake "$fish_target" \
    --snakefile workflow/fish_12s_host_depleted.smk --configfile "$fish_config" --cores "$SLURM_CPUS_PER_TASK" \
    --resources mem_mb=68000 --printshellcmds --rerun-incomplete --rerun-triggers mtime --latency-wait 60
rg --files -0 "$fish_root/union" | sort -z | xargs -0 sha256sum > "$run_record/union_output_sha256.txt"
