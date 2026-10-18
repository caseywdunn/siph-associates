#!/usr/bin/env bash
#SBATCH --job-name=coi-origins
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --output=logs/crustacean_coi/verification-%j.log
set -euo pipefail
module reset
run_record=data/results/crustacean_coi/executions/${SLURM_JOB_ID}
mkdir -p "$run_record"
exec > >(tee "$run_record/execution.log") 2>&1
trap 'code=$?; date -u > "$run_record/finished.txt"; printf "%s\n" "$code" > "$run_record/exit_code.txt"' EXIT
date -u > "$run_record/started.txt"
git rev-parse HEAD > "$run_record/git_head.txt"
git status --short > "$run_record/git_status.txt"
tar -czf "$run_record/source.tar.gz" scripts/crustacean_coi workflow/crustacean_coi*.smk workflow/rules/crustacean_coi*.smk config/crustacean_coi*.json dev_docs/crustacean-coi/validation/reference_review.tsv
sha256sum "$run_record/source.tar.gz" > "$run_record/source.sha256"
sha256sum data/results/crustacean_coi/panel/2026-10-17/{reference_metadata.tsv,target.fna,target.faa,annotated.faa} \
    data/results/crustacean_coi/reference/2026-10-17/mitochondrion.1.genomic.gbff.gz \
    data/results/crustacean_coi/host_depleted/2026-10-17/*/{candidates.fna,candidate_pairs.tsv,competitive.tsv} \
    data/results/crustacean_coi/validation/2026-10-17/heldout/{candidates.fna,candidate_pairs.tsv,competitive.tsv,truth.json} \
    data/results/crustacean_coi/validation/2026-10-17/panel/reference_metadata.tsv \
    data/results/crustacean_coi/host_depleted/2026-10-17/controls/truth.json \
    data/results/crustacean_coi/verification/2026-10-18/source/* > "$run_record/input_sha256.txt"
set -x
/vast/palmer/apps/avx2/software/DIAMOND/2.1.24-GCC-13.3.0/bin/diamond version
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python -c 'import Bio,sys; print(sys.version); print("Biopython",Bio.__version__)'
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python -m snakemake verification \
    --snakefile workflow/crustacean_coi_verification.smk --cores "$SLURM_CPUS_PER_TASK" \
    --resources mem_mb=18000 --printshellcmds --rerun-incomplete --rerun-triggers mtime --latency-wait 60
cp data/results/crustacean_coi/verification/2026-10-18/{output_sha256.txt,summary.json,gate.json} "$run_record/"
rg --files -0 data/results/crustacean_coi/verification/2026-10-18/{source,anchors,curated,queries,verification} | sort -z | xargs -0 sha256sum > "$run_record/verification_output_sha256.txt"
