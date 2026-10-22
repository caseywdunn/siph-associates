# Organisms associated with siphonophores

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23441460.svg)](https://doi.org/10.5281/zenodo.23441460)

This repository contains the analyses accompanying:

**Organisms associated with siphonophores revealed by short-read whole-genome sequencing**

*Manuscript citation and public link to be added.*

We examine DNA from bacteria, viruses, prey and parasites captured
alongside siphonophore DNA during genome sequencing. The study uses 213 specimens:
151 *Physalia*, 22 *Nanomia* and 40 other siphonophores. Start with the tables
below to connect the paper to its methods and supporting data;
understanding the analyses does not require running the code.

## Source studies

Each specimen is represented by one sequencing library. Counts below are the
unique libraries assigned to each primary study in this analysis; samples reused
across publications are counted once. The [library provenance table](data/metadata/library_provenance.tsv)
records shared samples among the Ahuja and Church studies; the
[Manko2026 specimen table](data/metadata/manko2026/sample_metadata.tsv) records
specimens from Mańko et al. (2026).

| Study | Samples in this analysis | Data | Role |
| --- | ---: | --- | --- |
| Ahuja et al. 2024, *Genome Biology and Evolution* ([paper](https://doi.org/10.1093/gbe/evae048)) | 33 | Genome skimming across siphonophore species | Broad sampling across the siphonophore phylogeny |
| Church et al. 2025, *Current Biology* ([paper](https://doi.org/10.1016/j.cub.2025.05.066)) | 151 | *Physalia* population samples and a chromosome-scale *P. physalis* genome | Intensive *Physalia* sampling and a host reference |
| Ahuja et al. 2026, *PLOS One* ([paper](https://doi.org/10.1371/journal.pone.0351247)) | 21 | *Nanomia* samples and a chromosome-scale *N. septata* genome | Intensive *Nanomia* sampling and a host reference |
| Mańko et al. 2026, *Zoological Journal of the Linnean Society* ([paper](https://doi.org/10.1093/zoolinnean/zlag150)) | 8 | *Bargmannia* genome-skimming libraries from Table S2 | Sampling across *Bargmannia* species |
| **Total** | **213** | | |

The [sample metadata](figures/cohort/sample_metadata.tsv) provide host names,
collection locations, tissues, depths and sequencing accessions.
The [host-depletion table](figures/host_depletion/mapping_stats.md)
([TSV](figures/host_depletion/mapping_stats.tsv)) reports reads ingested,
percent removed, and reads remaining for every library with a host reference.

## Finding analyses from the paper

The rows follow the current manuscript's Results headings in order. Figure and
table numbers refer to the manuscript and its separate supplementary materials.
Each guide explains the methods and links to their implementation; the data links
point to small tables retained in this repository.

| Results section in the paper | Figures and tables | Analyses and files |
| --- | --- | --- |
| Genome sequencing reveals a broad assemblage of associated organisms | Figure 1; Figures S1–S2; Table 1 | [Workflow overview](docs/workflow/README.md), [sample metadata](figures/cohort/sample_metadata.tsv), [sampling](figures/cohort/host_sampling.tsv), and [evidence summary](figures/overview/domain_summary.tsv) |
| Bacterial associates include diverse proteobacteria and recurrent Metamycoplasmataceae | Figure 2; Figure S5 | [Genome recovery](docs/workflow/phase3_assembly/README.md), [genome selection](docs/workflow/phase3_catalog/README.md), [detection criteria](docs/workflow/phase6_analysis/README.md), and [host distributions](figures/overview/microbial_host_distribution.tsv) |
| *Vibrio* diversity extends across *Physalia* and *Nanomia* | Text results | [Read mapping](docs/workflow/phase5_mapping/README.md), [assembly support](docs/workflow/phase6_analysis/README.md), and [bacterial detections](figures/cohort/bacterial_detections.tsv) |
| Metamycoplasmataceae have relatives in other marine animals | Figure 3; Figures S3 and S5 | [Genome phylogeny](docs/workflow/mycoplasmatales_phylogeny/README.md), [cnidarian marker comparisons](docs/workflow/cnidarian_16s/README.md), [Cassiopea comparisons](docs/workflow/cassiopea_16s/README.md), [tree data](figures/mycoplasmatales/source_phylogeny_focus.tsv), and [host incidence](figures/cohort/mycoplasmatales_incidence.tsv) |
| Metamycoplasmataceae differ in predicted ways of obtaining energy | Figure S6; Table S2 | [Gene functions](docs/workflow/mycoplasmatales_gene_content/README.md), [figure source data](figures/mycoplasmatales/source_gene_functions.tsv), and [nutritional annotations](figures/mdb_revision/amino_acid_modules.tsv) |
| Three bacterial lineages are frequent in southwestern Pacific and Indian Ocean samples | Figures S4 and S7; Table S1 | [Regional tests](docs/workflow/phase6_analysis/README.md), [matched host and regional comparisons](docs/coauthor_metadata_analysis.md), [incidence counts](figures/coauthor_metadata/targeted_incidence.tsv), and [C. acnes omission sensitivity](docs/cutibacterium_sensitivity.md) |
| Tissue and depth patterns reflect the sampling of colonies | Figures S3–S4 and S8 | [Tissue and depth methods](docs/coauthor_metadata_analysis.md), [tissue incidence](figures/coauthor_metadata/tissue_incidence.tsv), [depth design](figures/coauthor_metadata/host_depth_design.tsv), and [fish detections](figures/fish_12s/sample_summary.tsv) |
| Viral diversity includes viruses of associated bacteria | Table 1 | [Viral recovery and selection](docs/workflow/phase3_catalog/README.md), [predicted microbial hosts](docs/workflow/phase4_catalog/README.md), and [viral detections](figures/cohort/viral_detections.tsv) |
| Eukaryotic associates include prey and parasite lineages | Figure 4; Table 1 | [Ribosomal-marker analysis](docs/workflow/eukaryote_gate/README.md), [eukaryote evidence](figures/cohort/eukaryote_evidence.tsv), [fish 12S results](docs/fish_12s_results.md), [trematode comparison](docs/waki_sequence_comparison.md), and [Perkinsela bin review](docs/mag_domain_audit.md) |
| Initial comparisons recovered complementary evidence and some spurious assignments | Text and supplementary results | [Screening methods](docs/workflow/screens/README.md) and [screen comparisons](figures/screening/) |
| Siphonophore reference genomes greatly reduced the reads retained for associate analysis | Text and supplementary results | [Host-reference comparisons](docs/coauthor_metadata_analysis.md#host-reference-context), [current removal statistics](figures/coauthor_metadata/host_depletion_summary.tsv), and [historical fixed-catalog comparison](figures/coauthor_metadata/historical_host_concordance.tsv) |
| Human DNA and bacteria associated with handling | Text and supplementary results | [Human-read methods](docs/human_read_screen.md), [per-library screen](figures/human_screen/sample_summary.tsv), and [C. acnes comparison](docs/human_cutibacterium.md) |

[Co-detection summaries](figures/mdb_revision/cooccurrence.tsv) support Table S3
and the Discussion. Table S4 in the supplement inventories the accompanying data.
Figure 1 is an explanatory methods drawing maintained with the manuscript;
the executable workflow's generated rule graphs are linked below.

## Following the analysis steps

Read screening and assembly provide parallel evidence: every library is screened
and compared with shared microbial catalogs, even if it produces no assembled
microbial genome. Both host-depleted and other libraries contribute assemblies.
Host-read removal uses the *Physalia physalis* reference for all 151 *Physalia*
libraries, the *N. septata* reference for all ten *N. septata* libraries,
and the *Halistemma rubrum* reference for two *H. rubrum* libraries.
The other 50 libraries use no host reference.

| Analysis step | What it does | Methods and rule graphs |
| --- | --- | --- |
| Prepare reads | Quality-filter sequences; choose host-depleted reads or a smaller read subset for assembly | [Read preparation](docs/workflow/preprocess/README.md), [assembly inputs](docs/workflow/phase3_cohort/README.md), [graph](docs/workflow/graphs/inputs.svg) |
| Screen unassembled reads | Find candidate organisms and ribosomal markers before host-read removal | [Screening guide](docs/workflow/screens/README.md), [graph](docs/workflow/graphs/screening.svg) |
| Assemble and recover genomes | Reconstruct longer DNA sequences, microbial genomes and viral candidates from individual libraries | [Assembly guide](docs/workflow/phase3_assembly/README.md), [graph](docs/workflow/graphs/assembly.svg) |
| Build shared catalogs | Select recovered sequences, group close relatives and add nominated public microbial genomes | [Recovered sequences](docs/workflow/phase3_catalog/README.md), [shared catalogs](docs/workflow/phase4_catalog/README.md), [graph](docs/workflow/graphs/mapping_catalog.svg) |
| Detect associates across libraries | Map unassembled reads to the catalogs and assess support from the same specimen's assemblies | [Current mapping](docs/cohort_remapping.md), [detection criteria](docs/workflow/phase6_analysis/README.md), [mapping graph](docs/workflow/graphs/cohort_remapping.svg) |
| Identify other eukaryotes | Compare ribosomal-marker reads and assembled markers with reference sequences | [Eukaryote guide](docs/workflow/eukaryote_gate/README.md), [graph](docs/workflow/graphs/expanded_eukaryotes.svg), [trematode comparison](docs/waki_sequence_comparison.md) |
| Screen for fish prey | Search fish 12S in all 213 libraries: all retained host-depleted pairs for 151 *Physalia*, ten *N. septata* and two *H. rubrum* libraries; up to 200 million raw pairs before quality filtering for the other 50 | [Fish-screen methods and execution](docs/prey_marker_screen.md), [results](docs/fish_12s_results.md), [per-library read sources and counts](figures/fish_12s/sample_summary.tsv) |
| Audit human-compatible reads | Screen bounded raw-read prefixes with stringent paired alignments, simulated controls and candidate comparisons against broader references | [Methods and interpretation](docs/human_read_screen.md), [validation](dev_docs/human-screen/2026-10-22-validation.md) |
| Compare bacterial relationships and functions | Infer a genome tree, search earlier marker surveys and compare predicted gene functions | [Phylogeny](docs/workflow/mycoplasmatales_phylogeny/README.md), [tree graph](docs/workflow/graphs/phylogeny.svg), [gene functions](docs/workflow/mycoplasmatales_gene_content/README.md), [function graph](docs/workflow/graphs/gene_content.svg) |
| Compare distributions and assess sensitivity | Describe tissue and depth patterns, compare regions within hosts and hosts within regions, and assess sequencing-batch effects | [Statistical methods](docs/workflow/phase6_analysis/README.md), [graph](docs/workflow/graphs/expanded_statistics.svg), [collection context](docs/collection_context_analysis.md), [matched comparisons](docs/coauthor_metadata_analysis.md) |
| Produce figures and supporting tables | Summarize accepted results for the paper | [Current reporting guide](docs/expanded_reporting.md), [graph](docs/workflow/graphs/expanded_publication.svg), [figure and data index](figures/README.md) |

The [full workflow guide](docs/workflow/README.md) describes the analysis stages
and generated graphs. Each stage guide embeds its graph beside preview and run
commands. Graph arrows show computational dependencies; repeated sample jobs
are collapsed into one box per rule. Some reviewed inputs and standalone
analyses are documented separately, so the graphs do not represent every handoff.

## Data and interpretation

The [sample metadata](figures/cohort/sample_metadata.tsv) and
[detection tables](figures/cohort/) are starting points for exploring the results.
The [cohort selection](config/post_expansion_samples.tsv),
[input manifest](manifest.csv), [read inventory](data/metadata/raw_files.tsv)
and [source-study records](data/sources/) document the sequencing inputs.
Publication metadata retain updated host names alongside the names used in the
source studies. [Depth conventions](docs/metadata_depth_conventions.md) distinguish
recorded depths from assigned surface values.

All 213 libraries enter the descriptive summaries. Tests accounting for
sequencing batch use libraries sequenced on one physical flowcell; the generated
[flowcell table](figures/cohort/library_flowcells.tsv) records eligibility. The
Physalia geographic tests retain their 122 eligible libraries. Sampling, tissue, host reference availability,
geography and source study overlap. A detection establishes sequence evidence;
it does not by itself establish infection, symbiosis or host specificity.

## Reproducing an analysis

Start with [current cohort and reporting](docs/expanded_reporting.md), then the relevant stage guide and its
[environment and input instructions](docs/workflow/README.md#environment-and-inputs).
The workflow uses Snakemake to connect analysis steps and their files. Each guide
provides a separate **dry run** (a preview that executes no analysis jobs), a
launch command and a generated rule graph. Instructions for Yale's computing
cluster and job recovery are in the [operation guide](WORKFLOW.md).
The default Snakemake target checks input readiness; it does not run the whole paper.

This repository retains code, settings, small source files, figures and supporting
tables. Large sequencing files, reference databases and most sequence-analysis
outputs are stored separately; accessions and configured locations are described
in the input instructions. Some libraries lack public run accessions in the
recorded metadata, and some steps require reviewed intermediate inputs. A clone
alone is therefore insufficient for a complete rerun. The stage guides identify
these requirements and where outputs are written.

## Citation and development

Releases of this analysis repository are archived on Zenodo at
[doi:10.5281/zenodo.23441460](https://doi.org/10.5281/zenodo.23441460).
The manuscript citation and public supplementary-material link will be added when
available. [Workflow documentation checks](dev_docs/workflow-documentation.md)
and [execution records](RUNLOG.md) are maintained separately from the methods guides.
