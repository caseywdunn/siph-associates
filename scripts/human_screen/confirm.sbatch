#!/usr/bin/env bash
#SBATCH --job-name=human-candidate-nt
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/human-confirmation-%j.log
set -euo pipefail
module reset
module load BLAST+/2.15.0-gompi-2022b
human_scope=${1:-pilot}
human_record=data/results/human_screen/executions/${SLURM_JOB_ID}
mkdir -p "$human_record"
exec > >(tee "$human_record/execution.log") 2>&1
trap 'code=$?; date -u > "$human_record/finished.txt"; printf "%s\n" "$code" > "$human_record/exit_code.txt"' EXIT
date -u > "$human_record/started.txt"
git rev-parse HEAD > "$human_record/git_head.txt"
tar -czf "$human_record/source.tar.gz" scripts/human_screen workflow/human_confirmation.smk config/human_screen.json
sha256sum "$human_record/source.tar.gz" > "$human_record/source.sha256"
blastn -version
blastdbcmd -db /gpfs/gibbs/data/db/blast/nt -info > "$human_record/nt_database_info.txt"
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python -m snakemake all \
  --snakefile workflow/human_confirmation.smk --config scope="$human_scope" --cores "$SLURM_CPUS_PER_TASK" \
  --printshellcmds --rerun-incomplete --rerun-triggers mtime --latency-wait 60
