# Metamycoplasmataceae phylogeny figure

`phylogeny_focus.pdf` renders the 21-tip subtree containing all seven
siphonophore Metamycoplasmataceae representatives from the accepted
483-tip IQ-TREE 3.0.1 phylogeny (LG+F+I+R10). Plotting does not infer a new tree. MAGSP0012 remains in this
phylogenetic comparison despite exclusion of its source specimen from
biological incidence. MAGSP0031 is a divergent sister lineage of clade B.
The GTDB genus *Tyloplasma* is outside this figure's scope.

Build from the analysis repository root through SLURM:

```bash
sbatch -p ycga workflow/scripts/batch_mycoplasmatales_figures.sh
```

The batch launcher runs:

```bash
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python \
  workflow/scripts/plot_mycoplasmatales_figures.py
```

The existing Snakemake environment supplies Python, Biopython, NumPy and
Matplotlib. `manifest.json` records runtime versions, job ID, source revision,
the executed script and input/output SHA-256 checksums. Its source snapshot
identifies uncommitted plotting code. PDFs use embedded TrueType fonts;
PNG files are review previews.

## Sources and interpretation

The figure reads:

- `data/results/mycoplasmatales_phylogeny/iqtree3/2026-10-25/bac120.treefile`;
- `data/results/mycoplasmatales_phylogeny/iqtree3/2026-10-25/input/outgroup.txt`;
- `data/results/mycoplasmatales_phylogeny/tree_labels.tsv` for study and external labels;
- `config/mycoplasmatales_reference_taxonomy.tsv` for GTDB reference classifications.

The reference classifications retain the source identity recorded in the
[revision record](../../dev_docs/revision-2026-10-24-gene-story.md).
`source_phylogeny_focus.tsv` exports the displayed labels, source identities
and distances. The accepted [complete source tree](../../tables/mycoplasmatales/iqtree/bac120.treefile),
[IQ-TREE report](../../tables/mycoplasmatales/iqtree/bac120.iqtree),
[independent run trees](../../tables/mycoplasmatales/iqtree/bac120.runtrees),
[outgroup IDs](../../tables/mycoplasmatales/iqtree/outgroup.txt) and
[tree labels](../../tables/mycoplasmatales/tree_labels.tsv) are public
supporting files. The [input manifest](../../tables/mycoplasmatales/iqtree/inputs.json)
records sequence and source identities; the
[completion record](../../tables/mycoplasmatales/iqtree/inference_complete.json)
records inference integrity checks. The original
[FastTree output](../../tables/mycoplasmatales/mycoplasmatales.bac120.decorated.tree)
is retained as historical evidence and is not the source of the current figure.
Quality estimates in the original tree labels describe
phylogeny preparation; the plot does not display genome-quality values.

Node labels are SH-aLRT / ultrafast bootstrap support, both in percent,
from the best of five independent IQ-TREE runs. Black internal dots indicate
100% for both measures. `source_phylogeny_support.tsv` retains every displayed
node's exact support values, and `values.tex` supplies the reported values to
the manuscript. The weak clade-B/DT-68
relationship remains unresolved (64.2/71). The displayed subtree omits the
Acholeplasmataceae outgroup and does not establish relationships outside
the focal clade. Assembly GCA_929200685.1 appears as both an external
genome and a GTDB reference, explicitly marked `[dup.]`; these tips
represent the same genome.

## Export checks

Inspect the PNG preview and confirm embedded fonts with `pdffonts` before
copying the PDF and its source-table provenance into the manuscript.
The manuscript builds from its selected copies without sibling paths.
Earlier export checks are retained in the
[historical verification record](../../dev_docs/revision-2026-10-24-figure-history.md).
