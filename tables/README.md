# Supporting data for the siphonophore associates paper

These tables and selected marker sequences support the bacterial comparisons,
statistical tests and trematode comparisons described in the paper. Files are
copied unchanged from reviewed analysis outputs. The [export provenance](provenance.json)
records their original paths, sizes and SHA-256 checksums.

| Directory | Contents | Methods |
| --- | --- | --- |
| [mycoplasmatales](mycoplasmatales/) | Genome quality, tree labels, predicted functions, strict and relaxed ortholog/module assignments, and comparisons with earlier 16S surveys | [Phylogeny](../docs/workflow/mycoplasmatales_phylogeny/README.md), [gene content](../docs/workflow/mycoplasmatales_gene_content/README.md) |
| [phase6](phase6/) | Primary and sensitivity results, sequencing-flowcell eligibility, contamination tests, regional models, permutations, and phage-host links | [Statistical analyses](../docs/workflow/phase6_analysis/README.md) |
| [waki_comparison](waki_comparison/) | Trematode-marker alignments, competing reference matches, selected contigs and marker sequences, reference metadata, and search provenance | [Trematode comparisons](../docs/waki_sequence_comparison.md) |

In `phase6/`, `primary/` contains the primary analysis and `breadth_5pct/` and
`breadth_20pct/` contain detection-threshold sensitivity results. The
`sensitivity/` directory retains nomination-dependence and threshold summaries
alongside the historical fixed-catalog host-handling and read-capping controls.
Missing estimates and failed or untestable fits remain explicit.

In `mycoplasmatales/`, strict and relaxed assignments retain their original
thresholds. Incomplete genomes do not support claims of gene absence.
The `cnidarian_16s_` and `cassiopea_16s_` prefixes identify comparisons with
different surveys. Header-only match tables preserve an absence of qualifying
matches. The trematode files describe reference affinities, not new species
assignments. Selected FASTA files support those comparisons; they are not a
complete deposit of the reconstructed genomes or markers.

The phylogeny's `tree_labels.tsv` preserves genome-quality estimates recorded
when that tree was assembled. Use [genome_stats.tsv](mycoplasmatales/genome_stats.tsv)
for the quality values used in the manuscript and gene-function figure; join by
genome identifier rather than treating the tree-label estimates as current.
The Waki `reference_assembly_hits.tsv` likewise retains original host labels:
`Physalia sp1` for YPM-IZ-110931 is now *P. mikazuki*. Join specimen identifiers
to the [current sample metadata](../figures/cohort/sample_metadata.tsv) for current
host names. These retained source fields do not change the sequence comparisons.

For specimen metadata, detections and other figure source data, see the
[figure and data index](../figures/README.md).

## Verify or refresh the export

Python 3.10 or later and the standard library are sufficient. From the repository
root, verify a downloaded copy without the original analysis outputs:

```bash
python3 scripts/export_supporting_tables.py --check
```

To recreate the export where the original outputs are available, first preview
the reviewed [file selection](../config/supporting_tables.tsv):

```bash
python3 scripts/export_supporting_tables.py --dry-run
```

Then export and verify:

```bash
python3 scripts/export_supporting_tables.py
python3 scripts/export_supporting_tables.py --check
```

The exporter checks every source against the reviewed size and checksum before
copying. It does not rerun analyses. When reported results change, review and
update the selection before refreshing this export and the manuscript copies.
