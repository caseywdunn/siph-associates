#!/usr/bin/env bash
#SBATCH --job-name=phylogeny-hosts-build
#SBATCH --partition=ycga
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/phylogeny-hosts-build-%j.log
set -euo pipefail
module reset
module load poppler/22.12.0-GCC-12.2.0
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:$PATH"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
analysis_root="$PWD"
manuscript_root="/gpfs/ycga/work/dunn/cwd7/manuscript_siph_associates"
record="$PWD/data/results/manuscript_revision/2026-10-25/phylogeny_hosts/$SLURM_JOB_ID"
mkdir -p "$record/tmp" "$record/matplotlib"
export TMPDIR="$record/tmp"
export MPLCONFIGDIR="$record/matplotlib"
trap 'status=$?; printf "%s\n" "$status" > "$record/exit_code.txt"' EXIT
git rev-parse HEAD > "$record/git_head.txt"
cp workflow/scripts/plot_mycoplasmatales_figures.py \
  config/mycoplasmatales_source_hosts.tsv scripts/revision_phylogeny_hosts.sbatch "$record/"
cp figures/mycoplasmatales/source_phylogeny_focus.tsv "$record/phylogeny_before.tsv"
python workflow/scripts/plot_mycoplasmatales_figures.py
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
    for key in ('source', 'classification', 'duplicate_assembly',
                'distance_from_subtree_root', 'display_y'):
        assert before[tip][key] == after[tip][key], (tip, key)
    assert after[tip]['source_host'] != 'Not recorded'
    assert after[tip]['host_clade']
    assert after[tip]['source_host_record']
    assert after[tip]['source_host_url']
    if tip.startswith('MAGSP') or tip.startswith('EXT_') or tip == 'GB_GCA_929200685.1':
        assert before[tip]['source_host'] == after[tip]['source_host'], tip
for key in ('source_host', 'host_clade', 'source_host_record'):
    assert after['EXT_GCA_929200685.1'][key] == after['GB_GCA_929200685.1'][key]
print('All 21 tips, classifications, duplicate flags, branch distances and row positions unchanged.')
print('Existing host assignments preserved; every row has reviewed source and clade provenance.')
PY
pdffonts figures/mycoplasmatales/phylogeny_focus.pdf > "$record/figure3-fonts.txt"
cd "$manuscript_root"
python scripts/sync_assets.py --analysis-root "$analysis_root"
bash tmp/revision-commit-checks.sbatch
