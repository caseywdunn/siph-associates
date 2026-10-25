#!/usr/bin/env bash
#SBATCH --job-name=mycoplasmatales-iqtree
#SBATCH --partition=ycga
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/revision-iqtree-%j.log
set -euo pipefail
module reset
export PATH="/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:$PATH"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:?Submit this script through sbatch}"

run_dir="$PWD/data/results/mycoplasmatales_phylogeny/iqtree3/2026-10-25"
record="$run_dir/executions/${SLURM_JOB_ID:?}"
mkdir -p "$record" "$run_dir/tmp"
export TMPDIR="$run_dir/tmp"
trap 'status=$?; printf "%s\n" "$status" > "$record/exit_code.txt"' EXIT
date -u +%FT%TZ > "$record/started_utc.txt"
hostname > "$record/hostname.txt"
git rev-parse HEAD > "$record/git_head.txt"
git diff > "$record/source.diff"
cp scripts/revision_iqtree.sbatch workflow/scripts/prepare_mycoplasmatales_iqtree.py \
  envs/iqtree.yaml "$record/"
iqtree3 --version > "$record/iqtree_version.txt"
grep -q 'IQ-TREE version 3.0.1 ' "$record/iqtree_version.txt"
python --version > "$record/python_version.txt"
cp /gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/conda-meta/iqtree-3.0.1-h503566f_0.json \
  "$record/iqtree_conda_package.json"
sha256sum "$(command -v iqtree3)" > "$record/iqtree_binary.sha256"

# The original alignment and FastTree are immutable comparison inputs.
if [[ ! -e "$run_dir/input/inputs.json" ]]; then
  python workflow/scripts/prepare_mycoplasmatales_iqtree.py \
    --alignment data/results/mycoplasmatales_phylogeny/mycoplasmatales.bac120.msa.fasta.gz \
    --tree data/results/mycoplasmatales_phylogeny/mycoplasmatales.bac120.decorated.tree \
    --labels data/results/mycoplasmatales_phylogeny/tree_labels.tsv \
    --output-dir "$run_dir/input"
fi

(cd "$run_dir/input" && sha256sum --check source.sha256 generated.sha256)

# Reusing a checkpoint requires the same inputs and inference implementation.
if [[ -e "$run_dir/resume.sha256" ]]; then
  sha256sum --check "$run_dir/resume.sha256"
else
  sha256sum scripts/revision_iqtree.sbatch workflow/scripts/prepare_mycoplasmatales_iqtree.py \
    envs/iqtree.yaml \
    data/results/mycoplasmatales_phylogeny/mycoplasmatales.bac120.msa.fasta.gz \
    data/results/mycoplasmatales_phylogeny/mycoplasmatales.bac120.decorated.tree \
    data/results/mycoplasmatales_phylogeny/tree_labels.tsv \
    "$run_dir/input/alignment.fasta" "$run_dir/input/outgroup.txt" \
    "$(command -v iqtree3)" > "$run_dir/resume.sha256"
fi

outgroup=$(cat "$run_dir/input/outgroup.txt")
test -n "$outgroup"
set -x
iqtree3 -s "$run_dir/input/alignment.fasta" --seqtype AA \
  -m MFP --runs 5 -B 1000 --alrt 1000 --bnni \
  -T AUTO --threads-max "$SLURM_CPUS_PER_TASK" --mem 56G \
  --seed 20261025 -o "$outgroup" --prefix "$run_dir/bac120"
set +x

# Completion of inference is separate from reviewing topology and support.
python - "$run_dir" <<'PY'
import json
import sys
from pathlib import Path
from Bio import Phylo, SeqIO

root = Path(sys.argv[1])
expected = {r.id for r in SeqIO.parse(root / 'input/alignment.fasta', 'fasta')}
tree = Phylo.read(root / 'bac120.treefile', 'newick')
observed = [tip.name for tip in tree.get_terminals()]
if len(observed) != len(set(observed)) or set(observed) != expected:
    raise ValueError('IQ-TREE output tips do not match the complete input alignment')
run_trees = list(Phylo.parse(root / 'bac120.runtrees', 'newick'))
if len(run_trees) != 5:
    raise ValueError(f'Expected five independent run trees, found {len(run_trees)}')
for suffix in ('iqtree', 'log', 'treefile', 'runtrees'):
    path = root / f'bac120.{suffix}'
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'Missing IQ-TREE output: {path}')
(root / 'inference_complete.json').write_text(json.dumps({
    'tips': len(observed), 'independent_searches': len(run_trees),
    'topology_and_support_review': 'pending',
}, indent=2) + '\n')
PY
date -u +%FT%TZ > "$record/completed_utc.txt"
