# Figures and supporting data for the paper

These figures and small data tables summarize the organisms detected with
siphonophores and the relationships and predicted functions of their bacteria.
The [paper-to-analysis index](../README.md#finding-analyses-from-the-paper)
connects them to the Results sections. Figure numbers below refer to the current
manuscript and S1 Appendix. Supplementary figures use letters A–H; tables use A–D.

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
The retained `overview/study_overview.pdf` concept sketch has photo placeholders
and is not that manuscript figure.
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
| Fig A in S1 Appendix: host sampling | [Samples by host](cohort/sampling.pdf) | [Sampling counts](cohort/host_sampling.tsv), [sample metadata](cohort/sample_metadata.tsv) |
| Fig B in S1 Appendix: collection depth | [Depth and host context](collection_context/depth_host_overview.pdf) | [Specimen records](collection_context/library_context.tsv), [depth summaries](collection_context/depth_summary.tsv) |
| Fig C in S1 Appendix: occurrence of the DT-68 bacterial lineage | [Nanomia and Resomia comparison](collection_context/dt68_depth_context.pdf) | [Specimen records](collection_context/dt68_specimens.tsv), [depth summaries](collection_context/dt68_depth_summary.tsv) |
| Fig D in S1 Appendix: tissue and region | [Tissue and regional context](collection_context/tissue_context.pdf) | [Tissue counts](collection_context/tissue_summary.tsv), [regional counts](collection_context/tissue_region_incidence.tsv) |
| Fig E in S1 Appendix: recurrent mycoplasmas | [Mycoplasmatales incidence](cohort/mycoplasmatales_incidence.pdf) | [Positive and total specimens](cohort/mycoplasmatales_incidence.tsv) |
| Fig F in S1 Appendix: predicted functions | [Gene functions](mycoplasmatales/gene_functions.pdf) | [Function assignments](mycoplasmatales/source_gene_functions.tsv) |
| Fig G in S1 Appendix: matched host and region comparisons | [Targeted bacterial incidence](coauthor_metadata/targeted_incidence.pdf) | [Incidence counts](coauthor_metadata/targeted_incidence.tsv), [community tests](coauthor_metadata/community_tests.tsv), [tested specimens](coauthor_metadata/community_membership.tsv) |
| Fig H in S1 Appendix: depth and associate detections | [Depth patterns](coauthor_metadata/depth_patterns.pdf) | [Specimen context](coauthor_metadata/specimen_context.tsv), [depth summaries](coauthor_metadata/depth_incidence.tsv), [depth design](coauthor_metadata/host_depth_design.tsv) |

PDFs are vector figures; PNG files alongside them are previews. The
[cohort tables](cohort/) retain all 213 libraries, including those with no
accepted detection. The microbial-distribution cells count each specimen once
per displayed group. The eukaryote figure shows only specimens with accepted detections;
its broad taxonomic groups are not comparable measures of species richness.
The historically named `bacterial_detections.tsv` includes both bacterial and
archaeal catalog targets. `summary.json` reports their combined `microbial_*`
counts and separate `bacterial_*` and `archaeal_*` counts; the domain table also
separates them.
[Collection-context methods](../docs/collection_context_analysis.md) explain the
recorded versus assigned depths and the overlap of host, tissue and geography.

Table 1 uses the [domain summary](overview/domain_summary.tsv).
[Selected nutritional and defense annotations](../docs/mdb_revision.md) and
[co-detection counts](mdb_revision/cooccurrence.tsv) support Tables B and C in S1 Appendix.
The other supplementary tables summarize detection-threshold sensitivity
(Table A in S1 Appendix) and inventory the accompanying data (Table D); their manuscript
formatting is maintained with the supplement. See the
[statistical methods](../docs/workflow/phase6_analysis/README.md) for sensitivity
outputs and the [trematode comparison](../docs/waki_sequence_comparison.md) for
additional marker evidence.
The [supporting-data export](../tables/README.md) includes the complete selected
genome-quality, gene-content, statistical and trematode-comparison files.

The [metadata comparison bundle](coauthor_metadata/) also contains tissue
incidence, current host-depletion summaries and the historical fixed-catalog
read-handling comparison. [Methods and interpretation](../docs/coauthor_metadata_analysis.md)
describe which host and regional contrasts can be tested and where tissue and
depth remain confounded with sampling.

The [human-screen bundle](human_screen/) contains per-library prefix counts,
simulated controls and representative-sequence comparisons against broader
references. The [human/C. acnes bundle](human_cutibacterium/) gives the descriptive
specimen join, and the [C. acnes omission bundle](cutibacterium_sensitivity/)
reports the two community contrasts without that bacterial target. Each includes
generated manuscript values and provenance; the
[human-read guide](../docs/human_read_screen.md) explains the selected-prefix
denominators and identification limits.

## Rebuild figures from accepted results

Use [current cohort and downstream reporting](../docs/expanded_reporting.md)
for the selected inputs and complete reporting sequence. The
[reviewed reporting workflow](../docs/viral_reporting_review.md) rebuilds cohort,
collection-context and overview figures from completed evidence, while the
[metadata workflow](../docs/coauthor_metadata_analysis.md) produces the matched
host/region and depth panels. The retained bacterial tree and function panels
have their own [figure guide](mycoplasmatales/README.md).
These reporting stages require completed upstream results. They do not run
sequence assembly or read mapping.

The figure producers are:

- [build_cohort.py](build_cohort.py): sample metadata, sampling, detections and
  Mycoplasmatales incidence. Publication host names are applied to the frozen
  analysis metadata while the original study names are retained.
- [plot_mycoplasmatales_figures.py](../workflow/scripts/plot_mycoplasmatales_figures.py):
  the focused bacterial tree and predicted functions.
- [analyze_collection_context.py](../scripts/analyze_collection_context.py) and
  [plot_collection_context.py](plot_collection_context.py): depth, tissue and
  regional summaries and figures.
- [Metadata workflow](../workflow/coauthor_metadata.smk): Figs G–H in S1 Appendix,
  matched community comparisons and host-reference summaries, using completed
  evidence tables as described in the [metadata guide](../docs/coauthor_metadata_analysis.md).

The small human-screen and C. acnes reports have separate bounded workflows;
their [workflow index](../docs/workflow/README.md) links the methods and entry points.

### Biological overview figures

Figures 2 and 4 are generated by
[plot_associate_overview.py](plot_associate_overview.py), connected through
`expanded_overview` in the expanded reporting workflow. The script can also be
run separately from completed tables. First generate the cohort and collection-context tables with the
reviewed reporting workflow above. Then, from the analysis repository root, use the
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
and outputs. The metadata and human-related bundles have their own provenance:
[metadata comparisons](coauthor_metadata/provenance.json),
[human screen](human_screen/provenance.json),
[human/C. acnes comparison](human_cutibacterium/provenance.json), and
[C. acnes omission](cutibacterium_sensitivity/provenance.json).
The [supporting-data provenance](../tables/provenance.json) identifies the
original analysis paths and checksums for the additional tables and selected
marker sequences.

The repository's [graph-freshness CI workflow](../.github/workflows/rulegraphs.yml)
checks the generated workflow graphs in its configured data-access environment.
[Graph regeneration and checks](../docs/workflow/README.md#regenerate-and-check-the-graphs)
use the same generator locally. These checks concern dependency diagrams;
figure and scientific-result validation are recorded in the producer manifests.
