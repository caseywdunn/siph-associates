#!/usr/bin/env bash
#SBATCH --job-name=host-depletion-report
#SBATCH --partition=ycga
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=logs/host-depletion-report-%j.log
set -euo pipefail
module reset
export XDG_CACHE_HOME="$PWD/.cache"
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python -m snakemake \
  --snakefile workflow/host_depletion_report.smk report_host_depletion \
  --cores 1 --nolock --printshellcmds \
  --config host_depletion_require_complete=true
