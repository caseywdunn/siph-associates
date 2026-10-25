# Supporting data for the siphonophore associates paper

These tables and selected marker sequences support the bacterial comparisons,
statistical tests and trematode comparisons described in the paper. Files are
copied unchanged from reviewed analysis outputs. The [export provenance](provenance.json)
records their original paths, sizes and SHA-256 checksums.

| Directory | Contents | Methods |
| --- | --- | --- |
| [mycoplasmatales](mycoplasmatales/) | Tree labels and comparisons with earlier 16S surveys | [Phylogeny](../docs/workflow/mycoplasmatales_phylogeny/README.md), [16S comparisons](../docs/workflow/cnidarian_16s/README.md) |
| [phase6](phase6/) | Primary and sensitivity results, sequencing-flowcell eligibility, contamination tests, regional models and permutations | [Statistical analyses](../docs/workflow/phase6_analysis/README.md) |
| [waki_comparison](waki_comparison/) | Trematode-marker alignments, competing reference matches, selected contigs and marker sequences, reference metadata, and search provenance | [Trematode comparisons](../docs/waki_sequence_comparison.md) |

In `phase6/`, `primary/` contains the primary analysis and `breadth_5pct/` and
`breadth_20pct/` contain detection-threshold sensitivity results. The
`sensitivity/` directory retains nomination-dependence and threshold summaries
alongside the historical fixed-catalog host-handling and read-capping controls.
Missing estimates and failed or untestable fits remain explicit.

In `mycoplasmatales/`, the `cnidarian_16s_` and `cassiopea_16s_` prefixes identify comparisons with
different surveys. Header-only match tables preserve an absence of qualifying
matches. The trematode files describe reference affinities, not new species
assignments. Selected FASTA files support those comparisons; they are not a
complete deposit of the reconstructed genomes or markers.

The phylogeny's `tree_labels.tsv` preserves genome-quality estimates recorded
when that tree was assembled. Use the [reviewed MAG catalog](../figures/manuscript_revision/catalog/mag_catalog.tsv)
for the quality values associated with the mapping catalog; join by
genome identifier rather than treating the tree-label estimates as current.
The Waki source alignments can retain original host labels; the specimen
summary uses the current name *P. mikazuki* for YPM-IZ-110931. Join specimen identifiers
to the [current sample metadata](../figures/cohort/sample_metadata.tsv) for current
host names. These retained source fields do not change the sequence comparisons.

For specimen metadata, detections and other figure source data, see the
[figure and data index](../figures/README.md).

Biological detection summaries use the 209 retained specimens and omit
contaminant taxa. Full evidence tables retain flagged observations from
all 213 screened libraries. [Exclusion decisions](../config/excluded_libraries.tsv)
and the [contaminant list](../config/contaminant_taxa.tsv) record this scope.
For eukaryotes, assembly support includes reads + assembly and assembly-only
detections; the two counts are distinct.

The organism-name crosswalk remains pending. It will record proposed NCBI
names separately from the GTDB classifications used in the paper.
Genome deposit preparation retains edited-sequence quality estimates
separately from the original mapping catalog.

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
