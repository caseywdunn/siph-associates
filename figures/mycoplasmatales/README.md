# Mycoplasmatales publication figures

These figures render existing analysis outputs; they do not rerun annotations
or infer a tree. Build from the analysis repository root on a compute node:

```bash
sbatch workflow/scripts/batch_mycoplasmatales_figures.sh
```

The batch command is:

```bash
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python \
  workflow/scripts/plot_mycoplasmatales_figures.py
```

The existing Snakemake conda environment supplies Python, Biopython, NumPy and
Matplotlib. Exact runtime versions, SLURM job ID, source revision, executed
script snapshot, and input/output SHA-256 checksums are recorded in
`manifest.json` when rendering succeeds. The source snapshot identifies
uncommitted plotting code as well as the repository revision. PDF output uses
embedded TrueType fonts; PNG files are review previews.

## Retained outputs and source evidence

- `gene_functions.pdf`: seven species-representative MAGs across selected
  conserved orthologs. Teal denotes a strict-threshold assignment, orange an
  assignment only under the primary relaxed rule, and grey no detected
  assignment. Counts greater than one are printed. The two low-quality MAGs
  are marked as presence-only observations. Family bars separately show
  strict and relaxed assignment frequencies among 130 near-complete GTDB
  representatives. Sources:
  `data/results/mycoplasmatales_gene_content/summary/focal_functions.tsv` and
  `genome_stats.tsv` in the same directory. The plotted values are exported
  as `source_gene_functions.tsv`.
- `phylogeny_focus.pdf`: the unchanged 21-tip MRCA subtree containing all
  seven siphonophore MAGs from the retained 483-tip tree. Taxonomy and
  external host labels come from the retained tables, not inferred host
  associations. Sources:
  `data/results/mycoplasmatales_phylogeny/mycoplasmatales.bac120.decorated.tree`,
  `tree_labels.tsv` in that directory, and gene-content `genome_stats.tsv`.
  Tip labels, source identities, and plotted distances are exported as
  `source_phylogeny_focus.tsv`.

The complete source tree remains in the analysis results. The displayed subtree
does not show the Acholeplasmataceae outgroup or justify relationships outside
the focal clade. FastTree local support is not bootstrap support; the weak
clade-B/DT-68 relationship must not be described as resolved. Assembly
GCA_929200685.1 appears twice in the source tree, as an external genome and as a
GTDB reference; both are retained and explicitly marked `[dup.]`. These tips
are not independent genomes.

## Interpretation and validation

The relaxed tier is the documented primary annotation rule, not a validation
of every predicted function. Its calibration and limits are in
`docs/mycoplasmatales_gene_content_plan.md`. Gene detection does not establish
enzyme activity or host phenotype. No detected assignment does not establish
biological absence, especially for MAGSP0029 and MAGSP0031. The E1 paralogs
are labelled 2-oxoacid dehydrogenase, preserving substrate uncertainty.

The script checks selected genome counts, unique KO/tier records, nesting of
strict within relaxed assignments, unique tree tips, and retention of all
seven study MAGs. Before manuscript export, inspect both PNG previews and
check the PDFs with `pdffonts`. Record those checks here after the SLURM run.

Copy selected PDFs plus their manifest and source-table provenance into the
manuscript repository. The manuscript must not depend on sibling paths at
LaTeX build time.

### Export verification, 2026-10-03

Both plots were rendered in compute allocation 11793267 using the bounded
publication workflow. Visual review found a clipped second footnote line in
the first gene-function export; tight output bounds corrected it and both
exports were rerendered. Final previews show complete footnotes, readable
labels, the explicit duplicate assembly and the weak internal support.
`pdffonts` verified embedded subset TrueType fonts. Plotted values, source
hashes and exact producing code are recorded in `manifest.json` and
`provenance/`. Selected PDFs and source tables were copied into the manuscript
repository and inspected again in the compiled manuscript.
