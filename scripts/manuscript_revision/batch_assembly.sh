#!/usr/bin/env bash
#SBATCH --job-name=revision-assembly-support
#SBATCH --partition=ycga
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/revision-assembly-%A_%a.log
set -euo pipefail
module reset
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:/vast/palmer/apps/avx2/software/miniconda/24.7.1/condabin:$PATH"
revision_root="data/results/manuscript_revision/2026-10-24"
export TMPDIR="$PWD/$revision_root/scratch"
export XDG_CACHE_HOME="$PWD/$revision_root/cache"
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME"
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" config/revision_assembly_samples.txt)
    [[ -n "$sample" ]]
    target="$revision_root/phase6_analysis/v3/assembly_support/$sample.paf"
else
    target="$revision_root/phase6_analysis/v3/assembly_support/bacterial_catalog.asm5.mmi"
fi
python -m snakemake --snakefile workflow/manuscript_revision.smk "$target" \
    --cores "$SLURM_CPUS_PER_TASK" --nolock --rerun-incomplete --printshellcmds \
    --allowed-rules index_catalog_minimap2 align_assembly_catalog_minimap2
