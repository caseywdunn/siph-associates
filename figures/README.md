# Figures and supporting data for the paper

These figures and small data tables summarize the organisms detected with
siphonophores and the relationships and predicted functions of their bacteria.
The [paper-to-analysis index](../README.md#finding-analyses-from-the-paper)
connects them to the Results sections. Figure numbers below refer to the current
manuscript and its separate supplement.

The [sample metadata guide](../docs/sample_metadata.md) defines museum, BioSample,
collecting-event, tissue-sample, and analysis identifiers. It also links the
remaining metadata gaps and describes how to update this table independently
of sequence analyses.

The [host-depletion mapping table](host_depletion/mapping_stats.md)
([TSV](host_depletion/mapping_stats.tsv)) gives input reads, percent removed,
and retained reads for each reference-bearing sample. See the
[reporting workflow](../docs/host_depletion_report.md) to rebuild it from the
existing mapping summaries.

Figure 1 is an explanatory methods drawing maintained with the manuscript.
It is separate from the generated [Snakemake rule graphs](../docs/workflow/README.md#read-the-graphs),
which describe computational dependencies. The manuscript builds from selected
copies of figures and tables; it does not read files from this repository at
compilation time.

## Selected figures and source tables

| Paper figure | Figure file | Supporting data |
| --- | --- | --- |
| Figure 2: bacteria across hosts | [Microbial distributions](overview/microbial_host_distribution.pdf) | [Host counts and denominators](overview/microbial_host_distribution.tsv), [genome assignments](overview/microbial_display_assignments.tsv) |
| Figure 3: bacterial relationships | [Metamycoplasmataceae tree](mycoplasmatales/phylogeny_focus.pdf) | [Tree data and labels](mycoplasmatales/source_phylogeny_focus.tsv) |
| Figure 4: eukaryotic associates | [Eukaryote evidence](overview/eukaryote_evidence.pdf) | [Displayed evidence](overview/eukaryote_display.tsv), [full evidence table](cohort/eukaryote_evidence.tsv) |
| Figure S1: host sampling | [Samples by host](cohort/sampling.pdf) | [Sampling counts](cohort/host_sampling.tsv), [sample metadata](cohort/sample_metadata.tsv) |
| Figure S2: collection depth | [Depth and host context](collection_context/depth_host_overview.pdf) | [Specimen records](collection_context/library_context.tsv), [depth summaries](collection_context/depth_summary.tsv) |
| Figure S3: DT-68 occurrence | [Nanomia and Resomia comparison](collection_context/dt68_depth_context.pdf) | [Specimen records](collection_context/dt68_specimens.tsv), [depth summaries](collection_context/dt68_depth_summary.tsv) |
| Figure S4: tissue and region | [Tissue and regional context](collection_context/tissue_context.pdf) | [Tissue counts](collection_context/tissue_summary.tsv), [regional counts](collection_context/tissue_region_incidence.tsv) |
| Figure S5: recurrent mycoplasmas | [Mycoplasmatales incidence](cohort/mycoplasmatales_incidence.pdf) | [Positive and total specimens](cohort/mycoplasmatales_incidence.tsv) |
| Figure S6: predicted functions | [Gene functions](mycoplasmatales/gene_functions.pdf) | [Function assignments](mycoplasmatales/source_gene_functions.tsv) |

PDFs are vector figures; PNG files alongside them are previews. The
[cohort tables](cohort/) retain all 205 libraries, including those with no
accepted detection. The microbial-distribution cells count each specimen once
per displayed group. The eukaryote figure shows only the 39 positive specimens;
its broad taxonomic groups are not comparable measures of species richness.
[Collection-context methods](../docs/collection_context_analysis.md) explain the
recorded versus assigned depths and the overlap of host, tissue and geography.

Table 1 uses the [domain summary](overview/domain_summary.tsv).
[Selected nutritional and defense annotations](../docs/mdb_revision.md) and
[co-detection counts](mdb_revision/cooccurrence.tsv) support Tables S2 and S3.
The other supplementary tables summarize detection-threshold sensitivity
(Table S1) and inventory the accompanying data (Table S4); their manuscript
formatting is maintained with the supplement. See the
[statistical methods](../docs/workflow/phase6_analysis/README.md) for sensitivity
outputs and the [trematode comparison](../docs/waki_sequence_comparison.md) for
additional marker evidence.

## Rebuild figures from accepted results

Use the [publication workflow guide](../docs/workflow/publication_figures/README.md)
for the canonical environment, separate dry-run and launch commands, and embedded
rule graph. Its four reporting rules produce cohort figures, the bacterial tree
and function panels, and collection-context tables and plots. It restricts
execution to those reporting rules and requires completed upstream results.
It does not run sequence assembly, mapping or statistical inference.

The figure producers are:

- [build_cohort.py](build_cohort.py): sample metadata, sampling, detections and
  Mycoplasmatales incidence. Publication host names are applied to the frozen
  analysis metadata while the original study names are retained.
- [plot_mycoplasmatales_figures.py](../workflow/scripts/plot_mycoplasmatales_figures.py):
  the focused bacterial tree and predicted functions.
- [analyze_collection_context.py](../scripts/analyze_collection_context.py) and
  [plot_collection_context.py](plot_collection_context.py): depth, tissue and
  regional summaries and figures.

### Biological overview figures

Figures 2 and 4 are generated by the standalone
[plot_associate_overview.py](plot_associate_overview.py), outside the reporting
rule graph. First generate the cohort and collection-context tables with the
publication workflow above. Then, from the analysis repository root, use the
existing Yale plotting environment:

```bash
MPLCONFIGDIR=/tmp/siph-associates-matplotlib \
  /gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python \
  figures/plot_associate_overview.py
```

This reads the accepted sample metadata, bacterial, viral and eukaryote detections,
and host summaries. It writes figures and source tables under `figures/overview/`.
The recorded plotting environment uses Python 3.11.4, Matplotlib 3.7.1 and NumPy
1.24.2; manifests record the runtime used for each export. Large sequence inputs
are not needed for these plots, but the accepted upstream tables are required.
There is no separate dry-run mode for this plotting script.

## Source identity and checks

[Cohort checksums](cohort/source_checksums.tsv), the
[tree/function manifest](mycoplasmatales/manifest.json),
[collection-table provenance](collection_context/provenance.json),
[collection-plot manifest](collection_context/plot_manifest.json), and
[overview manifest](overview/manifest.json) identify the inputs, producing code
and outputs. The manuscript's transfer record identifies the selected copied
assets and their checksums. This separates figure presentation from sequence
processing and preserves the source of reported results.

The repository's [graph-freshness CI workflow](../.github/workflows/rulegraphs.yml)
checks the generated workflow graphs in its configured data-access environment.
[Graph regeneration and checks](../docs/workflow/README.md#regenerate-and-check-the-graphs)
use the same generator locally. These checks concern dependency diagrams;
figure and scientific-result validation are recorded in the producer manifests.
