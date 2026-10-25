# Metamycoplasmataceae phylogeny figure

`phylogeny_focus.pdf` renders the 21-tip subtree containing all seven
siphonophore Metamycoplasmataceae representatives from the retained
483-tip phylogeny. It does not infer a new tree. MAGSP0012 remains in this
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

- `data/results/mycoplasmatales_phylogeny/mycoplasmatales.bac120.decorated.tree`;
- `data/results/mycoplasmatales_phylogeny/tree_labels.tsv` for study and external labels;
- `config/mycoplasmatales_reference_taxonomy.tsv` for GTDB reference classifications.

The reference classifications retain the source identity recorded in the
[revision record](../../dev_docs/revision-2026-10-24-gene-story.md).
`source_phylogeny_focus.tsv` exports the displayed labels, source identities
and distances. The [complete source tree](../../tables/mycoplasmatales/mycoplasmatales.bac120.decorated.tree)
and [tree labels](../../tables/mycoplasmatales/tree_labels.tsv) are public
supporting files. Quality estimates in the original tree labels describe
phylogeny preparation; the plot does not display genome-quality values.

Node labels are FastTree SH-like local supports. The weak clade-B/DT-68
relationship remains unresolved. The displayed subtree omits the
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
