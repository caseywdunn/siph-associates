#!/usr/bin/env bash
#SBATCH --job-name=fish12s-full-reads
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/fish_12s/full-reads-%A_%a.log
set -euo pipefail
module reset
module load BLAST+/2.15.0-gompi-2022b Bowtie2/2.5.1-GCC-12.2.0 SeqKit/2.8.1
export PATH=/gpfs/gibbs/project/dunn/cwd7/conda_envs/phyloflash_env/bin:$PATH
export PATH="$PATH:/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:/vast/palmer/apps/avx2/software/miniconda/24.7.1/condabin"
export XDG_CACHE_HOME="$PWD/.cache"
python=/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python
fish_config=${1:-config/fish_12s_full_reads.json}
task_index=${SLURM_ARRAY_TASK_ID:-0}
readarray -t selected < <("$python" - "$fish_config" "$task_index" <<'PY'
import csv, json, sys
from pathlib import Path
config = json.loads(Path(sys.argv[1]).read_text())
with Path(config['libraries']).open() as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
row = rows[int(sys.argv[2])]
print(row['analysis_sample_id'])
print(config['run'])
print(config['reference_metadata'])
print(config['libraries'])
PY
)
sample=${selected[0]}
fish_root=${selected[1]}
metadata=${selected[2]}
libraries=${selected[3]}
run_record="data/results/fish_12s/executions/${SLURM_JOB_ID}_${task_index}_${sample}"
mkdir -p "$run_record"
exec > >(tee "$run_record/execution.log") 2>&1
trap 'code=$?; date -u > "$run_record/finished.txt"; printf "%s\n" "$code" > "$run_record/exit_code.txt"' EXIT
date -u > "$run_record/started.txt"
git rev-parse HEAD > "$run_record/git_head.txt"
git status --short > "$run_record/git_status.txt"
tar -czf "$run_record/source.tar.gz" scripts/fish_12s workflow/fish_12s_full_reads.smk \
    workflow/rules/fish_12s_search.smk workflow/scripts/{trim_reads.py,common.py} \
    "$fish_config" "$libraries" envs/fastp.yaml
sha256sum "$run_record/source.tar.gz" "$metadata" > "$run_record/input_sha256.txt"
sha256sum -c data/results/fish_12s/positive_references/2026-10-17/output_sha256.txt
"$python" scripts/fish_12s/check_cohort_controls.py \
    --root data/results/fish_12s/cohort_validation/2026-10-18 \
    --truth data/results/fish_12s/positive_aligner_controls/2026-10-17/synthetic/known_answer/checks.json \
    --reference-metadata "$metadata" --output "$run_record/CONTROLS_PASSED"
blastn -version
bowtie2 --version
seqkit version
bbduk.sh --version
data/results/fish_12s/software/sortmerna_env/bin/sortmerna --version
set -x
# One disjoint library per array task. No original analysis rules are included.
"$python" -m snakemake "$fish_root/review/$sample/sample_summary.tsv" \
    --snakefile workflow/fish_12s_full_reads.smk --configfile "$fish_config" \
    --config config_path="$fish_config" --cores "$SLURM_CPUS_PER_TASK" --resources mem_mb=60000 \
    --use-conda --conda-frontend mamba --printshellcmds --rerun-incomplete \
    --rerun-triggers mtime --latency-wait 60 --nolock
rg --files -0 "$fish_root/union/$sample" "$fish_root/review/$sample" | sort -z | xargs -0 sha256sum \
    > "$run_record/output_sha256.txt"
