#!/usr/bin/env bash
#SBATCH --job-name=revision-figures-build
#SBATCH --partition=ycga
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/revision-figures-build-%j.log
set -euo pipefail
module reset
module load poppler/22.12.0-GCC-12.2.0
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:$PATH"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
analysis_root="$PWD"
manuscript_root="/gpfs/ycga/work/dunn/cwd7/manuscript_siph_associates"
record="$PWD/data/results/manuscript_revision/2026-10-25/figures/$SLURM_JOB_ID"
mkdir -p "$record/tmp" "$record/matplotlib"
export TMPDIR="$record/tmp"
export MPLCONFIGDIR="$record/matplotlib"
trap 'status=$?; printf "%s\n" "$status" > "$record/exit_code.txt"' EXIT
git rev-parse HEAD > "$record/git_head.txt"
cp figures/plot_associate_overview.py workflow/scripts/plot_mycoplasmatales_figures.py \
  scripts/revision_figures.sbatch "$record/"
sha256sum figures/overview/*.tsv > "$record/overview_before.sha256"
cp figures/mycoplasmatales/source_phylogeny_focus.tsv "$record/phylogeny_before.tsv"
python figures/plot_associate_overview.py
python workflow/scripts/plot_mycoplasmatales_figures.py
sha256sum --check "$record/overview_before.sha256" > "$record/overview_unchanged.txt"
python - "$record" <<'PY'
import csv
import sys
from pathlib import Path

record = Path(sys.argv[1])
def read(path):
    with path.open() as handle:
        return {row['tip']: row for row in csv.DictReader(handle, delimiter='\t')}
before = read(record / 'phylogeny_before.tsv')
after = read(Path('figures/mycoplasmatales/source_phylogeny_focus.tsv'))
assert before.keys() == after.keys()
for tip in before:
    for key in ('source', 'classification', 'duplicate_assembly', 'distance_from_subtree_root'):
        assert before[tip][key] == after[tip][key], (tip, key)
print('All tree tips, classifications, duplicate flags and branch distances unchanged.')
PY
pdffonts figures/mycoplasmatales/phylogeny_focus.pdf > "$record/figure3-fonts.txt"
pdffonts figures/overview/eukaryote_evidence.pdf > "$record/figure4-fonts.txt"
cd "$manuscript_root"
python scripts/sync_assets.py --analysis-root "$analysis_root"
bash tmp/revision-commit-checks.sbatch
