#!/usr/bin/env bash
#SBATCH --job-name=coi-validation
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --output=logs/crustacean_coi/validation-%j.log
set -euo pipefail
module reset
module load BLAST+/2.15.0-gompi-2022b
export PATH=/gpfs/gibbs/project/dunn/cwd7/conda_envs/phyloflash_env/bin:$PATH
run_record=data/results/crustacean_coi/executions/${SLURM_JOB_ID}
mkdir -p "$run_record"
exec > >(tee "$run_record/execution.log") 2>&1
trap 'code=$?; date -u > "$run_record/finished.txt"; printf "%s\n" "$code" > "$run_record/exit_code.txt"' EXIT
date -u > "$run_record/started.txt"
git rev-parse HEAD > "$run_record/git_head.txt"
git status --short > "$run_record/git_status.txt"
tar -czf "$run_record/source.tar.gz" scripts/crustacean_coi workflow/crustacean_coi*.smk workflow/rules/crustacean_coi*.smk config/crustacean_coi*.json config/crustacean_coi_pilot.tsv
sha256sum "$run_record/source.tar.gz" > "$run_record/source.sha256"
sha256sum data/results/crustacean_coi/panel/2026-10-17/{reference_metadata.tsv,target.fna,target.faa,competitive.fna} data/results/crustacean_coi/host_depleted/2026-10-17/*/{candidates.fna,candidate_pairs.tsv,competitive.tsv,recruitment.json} /gpfs/ycga/work/dunn/cwd7/databases/diamond/swissprot.dmnd > "$run_record/input_sha256.txt"
sha256sum -c data/results/fish_12s/positive_references/2026-10-17/output_sha256.txt
set -x
/vast/palmer/apps/avx2/software/DIAMOND/2.1.24-GCC-13.3.0/bin/diamond version
blastn -version
bbduk.sh --version
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python -m snakemake "${1:-validation}" \
    --snakefile workflow/crustacean_coi_validation.smk --cores "$SLURM_CPUS_PER_TASK" \
    --resources mem_mb=18000 --printshellcmds --rerun-incomplete --rerun-triggers mtime --latency-wait 60
cp data/results/crustacean_coi/validation/2026-10-17/{output_sha256.txt,sensitivity_summary.json,design.json} "$run_record/"

if test -s data/results/crustacean_coi/validation/2026-10-17/protein_competition/summary.json; then
    sha256sum data/results/crustacean_coi/validation/2026-10-17/protein_competition/* > "$run_record/protein_output_sha256.txt"
    cp data/results/crustacean_coi/validation/2026-10-17/protein_competition/summary.json "$run_record/protein_summary.json"
fi
