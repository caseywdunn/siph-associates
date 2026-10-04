# Publication figures

Figures are generated here from the analysis repository's accepted results.
Keep scripts, source tables and provenance together; copy selected small PDF
exports to the manuscript repository with its `scripts/sync_assets.py` helper.
The manuscript TeX build does not read from this repository.

## Builds

The bounded downstream correction/figure workflow can be submitted through
the existing persistent SLURM controller:

```bash
bash scripts/submit_review_corrections.sh
```

Inside an existing compute allocation (at least one CPU and 5 GB), run:

```bash
bash scripts/run_review_corrections.sh --dry-run
bash scripts/run_review_corrections.sh
```

These entry points whitelist downstream statistics, eukaryote reporting and
publication figure rules. Existing assembly and mapping products are inputs;
missing upstream files cause failure rather than sequence reanalysis.
On 2026-10-03, the workflow ran inside allocation 11793267 on r209u08n04,
using one CPU and a 4-GB workflow memory budget. Source base revision was
`0a7c3e85a0e8f4856a16e5fc75e6a0ec73ba7b19`, with review corrections preserved
as working-tree snapshots because Git metadata was read-only.

## Exports

| Export | Build script | Inputs and interpretation |
| --- | --- | --- |
| `cohort/sampling.pdf` | `build_cohort.py` | Deduplicated manifest; all 205 libraries and all 37 host labels |
| `cohort/mycoplasmatales_incidence.pdf` | `build_cohort.py` | Catalog v2 and unchanged primary grades; full descriptive denominators |
| `cohort/eukaryote_evidence.pdf` | `build_cohort.py` | Corrected eukaryote reporting; 61 nonredundant detections, not species richness |
| `mycoplasmatales/phylogeny_focus.pdf` | `../workflow/scripts/plot_mycoplasmatales_figures.py` | Existing bac120 tree and labels; branch lengths and local support preserved |
| `mycoplasmatales/gene_functions.pdf` | same | Existing KO tables; strict and relaxed-only assignments distinguished |

PDFs are vector exports with embedded TrueType fonts. PNG files are review
previews. `cohort/source_checksums.tsv` and `cohort/source_revision.txt` identify
the inputs and plotting code. `mycoplasmatales/manifest.json` also records
runtime package versions and output checksums; its `provenance/` directory
contains the exact plotting script. The manuscript copy manifest additionally
records hashes of every copied artifact and exact dirty-source snapshots.

All libraries, including those sequenced on multiple flowcells, contribute to
descriptive exports. Only flowcell-based inference uses the explicitly
identified single-flowcell subset. The separate figure directory README records
Mycoplasmatales interpretation and visual checks.
