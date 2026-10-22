#!/usr/bin/env bash
#SBATCH --job-name=human-read-audit
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --output=logs/human-screen-%j.log
set -euo pipefail
module reset
module load fastp/0.23.2-GCCcore-10.2.0
module load Bowtie2/2.5.1-GCC-12.2.0
export PATH=/gpfs/gibbs/project/dunn/cwd7/conda_envs/assembly_env/bin:$PATH
human_target=${1:-pilot}
human_memory_mb=$(( ${SLURM_MEM_PER_NODE:-40000} - 2000 ))
human_record=data/results/human_screen/executions/${SLURM_JOB_ID}
mkdir -p "$human_record"
exec > >(tee "$human_record/execution.log") 2>&1
trap 'code=$?; date -u > "$human_record/finished.txt"; printf "%s\n" "$code" > "$human_record/exit_code.txt"' EXIT
date -u > "$human_record/started.txt"
git rev-parse HEAD > "$human_record/git_head.txt"
git status --short > "$human_record/git_status.txt"
tar -czf "$human_record/source.tar.gz" scripts/human_screen workflow/human_screen.smk config/human_screen.json config/post_expansion_samples.tsv
sha256sum "$human_record/source.tar.gz" > "$human_record/source.sha256"
bowtie2 --version
samtools --version
fastp --version
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python scripts/human_screen/record_references.py \
  --config config/human_screen.json --output "$human_record/reference_identity.json"
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python -m snakemake "$human_target" \
  --snakefile workflow/human_screen.smk --cores "$SLURM_CPUS_PER_TASK" --resources mem_mb="$human_memory_mb" \
  --printshellcmds --rerun-incomplete --rerun-triggers mtime --latency-wait 60
