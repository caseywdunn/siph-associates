#!/usr/bin/env bash
# Explicit large-memory launch after reviewing preparation and adaptor reports.
#SBATCH --job-name=deposit-fcs-gx
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=720G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/deposit-fcs-gx-%j.log
set -euo pipefail
module reset
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:$PATH"
export TMPDIR="$PWD/data/results/deposit/tmp/fcsgx-controller-$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
record="$PWD/data/results/deposit/executions/fcsgx-$SLURM_JOB_ID"
mkdir -p "$record"
git rev-parse HEAD > "$record/git_head.txt"
git diff > "$record/source.diff"
cp config/deposit_preparation.json config/mag_organism_names.tsv config/fcs_gx_division_overrides.tsv workflow/rules/deposit_preparation.smk workflow/scripts/check_deposit_fcsgx_inputs.py workflow/scripts/package_deposit_sequences.py workflow/scripts/prepare_deposit_mag.py "$record/"
sha256sum "$record/"* > "$record/source.sha256"
trap 'status=$?; printf "%s\n" "$status" > "$record/exit_code.txt"' EXIT
# Memory resources permit one screen at a time, reusing the GPFS database in the
# node's page cache. The database is never copied to /tmp or /dev/shm.
python -m snakemake -s workflow/deposit_preparation.smk deposit_fcsgx \
  --cores 16 --resources mem_mb=720000 --config enable_fcsgx=true \
  --rerun-triggers mtime --rerun-incomplete --printshellcmds "$@"
