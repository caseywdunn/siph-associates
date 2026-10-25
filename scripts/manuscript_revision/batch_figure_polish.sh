#!/usr/bin/env bash
#SBATCH --job-name=revision-figure-polish
#SBATCH --partition=ycga
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/revision-figure-polish-%j.log
set -euo pipefail
module reset
module load poppler/22.12.0-GCC-12.2.0
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:$PATH"
record="data/results/manuscript_revision/2026-10-24/figure_polish/$SLURM_JOB_ID"
mkdir -p "$record/tmp" "$record/matplotlib" "$record/cache"
export TMPDIR="$PWD/$record/tmp"
export MPLCONFIGDIR="$PWD/$record/matplotlib"
export XDG_CACHE_HOME="$PWD/$record/cache"
exec > >(tee "$record/execution.log") 2>&1
trap 'code=$?; date -u > "$record/finished.txt"; printf "%s\n" "$code" > "$record/exit_code.txt"' EXIT
date -u > "$record/started.txt"
git rev-parse HEAD > "$record/git_head.txt"
git status --short > "$record/git_status.txt"
cp scripts/summarize_contamination_revision.py workflow/scripts/analyze_coauthor_communities.R "$record/"
cp "$0" "$record/batch_figure_polish.sh"
sha256sum "$record/"*.py "$record/"*.R "$record/"*.sh > "$record/source.sha256"
sha256sum \
  figures/coauthor_metadata/community_tests.tsv \
  figures/coauthor_metadata/dispersion_tests.tsv \
  figures/coauthor_metadata/community_membership.tsv \
  figures/coauthor_metadata/community_values.json \
  figures/coauthor_metadata/values.tex \
  figures/contamination/sample_screen.tsv \
  figures/contamination/bacterial_contaminant_detections.tsv \
  figures/contamination/eukaryote_contaminant_detections.tsv \
  figures/contamination/contaminant_targets.tsv \
  figures/contamination/summary.json > "$record/scientific_before.sha256"
python scripts/summarize_contamination_revision.py
revision_rscript="/gpfs/gibbs/project/dunn/cwd7/conda_envs/stats_env/bin/Rscript"
"$revision_rscript" -e 'stopifnot(capabilities("cairo"))'
"$revision_rscript" workflow/scripts/analyze_coauthor_communities.R \
  config/coauthor_metadata_analysis.json --plots-only
python scripts/record_coauthor_metadata.py --config config/coauthor_metadata_analysis.json
sha256sum --check "$record/scientific_before.sha256" > "$record/scientific_unchanged.txt"
python - "$record" <<'PY'
import json
import subprocess
import sys
from pathlib import Path

result = {}
for name in ('depth_patterns', 'targeted_incidence'):
    path = Path('figures/coauthor_metadata') / f'{name}.pdf'
    output = subprocess.check_output(['pdffonts', str(path)], text=True)
    Path(sys.argv[1], f'{name}.fonts.txt').write_text(output)
    lines = output.splitlines()
    column = lines[0].index('emb')
    fonts = [line for line in lines[2:] if line.strip()]
    if not fonts or any(line[column:column+3] != 'yes' for line in fonts):
        raise ValueError(f'Unembedded or missing fonts: {path}')
    result[str(path)] = {'embedded_fonts': len(fonts), 'status': 'PASS'}
sys.path.insert(0, '/gpfs/ycga/work/dunn/cwd7/manuscript_siph_associates/scripts')
import sync_assets
result['figure_provenance'] = {
    'verified_files': len(sync_assets.recorded_figure_hashes(Path.cwd())),
    'status': 'PASS',
}
Path(sys.argv[1], 'validation.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
PY
