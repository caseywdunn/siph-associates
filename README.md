# Organisms associated with siphonophores

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23441460.svg)](https://doi.org/10.5281/zenodo.23441460)

This repository contains the analyses accompanying:

**Organisms associated with siphonophores revealed by short-read whole-genome sequencing**

*Manuscript citation and public link to be added.*

We examine DNA from bacteria, archaea, viruses, prey and parasites captured
alongside siphonophore DNA during genome sequencing. The study uses 205 specimens:
151 *Physalia*, 22 *Nanomia* and 32 other siphonophores. Start with the tables
below to connect the paper to its methods and supporting data;
understanding the analyses does not require running the code.

## Source studies

Each specimen is represented by one sequencing library. Counts below are the
unique libraries assigned to each primary study in this analysis; samples reused
across publications are counted once. The [library provenance table](data/metadata/library_provenance.tsv)
retains their relationships to all source studies.

| Study | Samples in this analysis | Data | Role |
| --- | ---: | --- | --- |
| Ahuja et al. 2024, *Genome Biology and Evolution* ([paper](https://doi.org/10.1093/gbe/evae048)) | 33 | Genome skimming across siphonophore species | Broad sampling across the siphonophore phylogeny |
| Church et al. 2025, *Current Biology* ([paper](https://doi.org/10.1016/j.cub.2025.05.066)) | 151 | *Physalia* population samples and a chromosome-scale *P. physalis* genome | Intensive *Physalia* sampling and a host reference |
| Ahuja et al. 2026, *PLOS One* ([paper](https://doi.org/10.1371/journal.pone.0351247)) | 21 | *Nanomia* samples and a chromosome-scale *N. septata* genome | Intensive *Nanomia* sampling and a host reference |
| **Total** | **205** | | |

One *Nanomia* library comes from Ahuja et al. 2024, giving 22 across studies.
The [sample metadata](figures/cohort/sample_metadata.tsv) provide host names,
collection locations, tissues, depths and sequencing accessions.

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
| Three bacterial lineages are frequent in southwestern Pacific and Indian Ocean samples | Figure S4; Table S1 | [Regional tests and sensitivity analyses](docs/workflow/phase6_analysis/README.md), [depth and tissue comparisons](docs/collection_context_analysis.md), and [regional detection counts](figures/collection_context/tissue_region_incidence.tsv) |
| One archaeal lineage recurs in *Physalia* and *Nanomia* | Figure 2; Table 1 | [Detection methods](docs/workflow/phase6_analysis/README.md), [bacterial/archaeal detections](figures/cohort/bacterial_detections.tsv), and [host distributions](figures/overview/microbial_host_distribution.tsv) |
| Viral diversity includes viruses of associated bacteria | Table 1 | [Viral recovery and selection](docs/workflow/phase3_catalog/README.md), [predicted microbial hosts](docs/workflow/phase4_catalog/README.md), and [viral detections](figures/cohort/viral_detections.tsv) |
| Eukaryotic associates include prey and parasite lineages | Figure 4; Table 1 | [Ribosomal-marker analysis](docs/workflow/eukaryote_gate/README.md), [eukaryote evidence](figures/cohort/eukaryote_evidence.tsv), and [direct trematode comparison](docs/waki_sequence_comparison.md) |

[Co-detection summaries](figures/mdb_revision/cooccurrence.tsv) support Table S3
and the Discussion. Table S4 in the supplement inventories the accompanying data.
Figure 1 is an explanatory methods drawing maintained with the manuscript;
the executable workflow's generated rule graphs are linked below.

## Following the analysis steps

Read screening and assembly provide parallel evidence: every library is screened
and compared with shared microbial catalogs, even if it produces no assembled
microbial genome. Both host-depleted and other libraries contribute assemblies.
Host-read removal uses the *Physalia physalis* reference for all 151 *Physalia*
libraries and the *N. septata* reference for all ten *N. septata* libraries;
the other 44 use no host reference.

| Analysis step | What it does | Methods and rule graphs |
| --- | --- | --- |
| Prepare reads | Quality-filter sequences; choose host-depleted reads or a smaller read subset for assembly | [Read preparation](docs/workflow/preprocess/README.md), [assembly inputs](docs/workflow/phase3_cohort/README.md), [graph](docs/workflow/graphs/inputs.svg) |
| Screen unassembled reads | Find candidate organisms and ribosomal markers before host-read removal | [Screening guide](docs/workflow/screens/README.md), [graph](docs/workflow/graphs/screening.svg) |
| Assemble and recover genomes | Reconstruct longer DNA sequences, microbial genomes and viral candidates from individual libraries | [Assembly guide](docs/workflow/phase3_assembly/README.md), [graph](docs/workflow/graphs/assembly.svg) |
| Build shared catalogs | Select recovered sequences, group close relatives and add nominated public microbial genomes | [Recovered sequences](docs/workflow/phase3_catalog/README.md), [shared catalogs](docs/workflow/phase4_catalog/README.md), [graph](docs/workflow/graphs/mapping_catalog.svg) |
| Detect associates across libraries | Map unassembled reads to the catalogs and assess support from the same specimen's assemblies | [Mapping](docs/workflow/phase5_mapping/README.md), [detection criteria](docs/workflow/phase6_analysis/README.md), [mapping graph](docs/workflow/graphs/mapping.svg) |
| Identify other eukaryotes | Compare ribosomal-marker reads and assembled markers with reference sequences | [Eukaryote guide](docs/workflow/eukaryote_gate/README.md), [graph](docs/workflow/graphs/eukaryotes.svg), [trematode comparison](docs/waki_sequence_comparison.md) |
| Compare bacterial relationships and functions | Infer a genome tree, search earlier marker surveys and compare predicted gene functions | [Phylogeny](docs/workflow/mycoplasmatales_phylogeny/README.md), [tree graph](docs/workflow/graphs/phylogeny.svg), [gene functions](docs/workflow/mycoplasmatales_gene_content/README.md), [function graph](docs/workflow/graphs/gene_content.svg) |
| Compare distributions and assess sensitivity | Describe host and tissue patterns, test regional variation and assess sequencing-batch effects | [Statistical methods](docs/workflow/phase6_analysis/README.md), [graph](docs/workflow/graphs/analysis.svg), [collection context](docs/collection_context_analysis.md) |
| Produce figures and supporting tables | Summarize accepted results for the paper | [Reporting guide](docs/workflow/publication_figures/README.md), [graph](docs/workflow/graphs/publication.svg), [figure and data index](figures/README.md) |

The [full workflow guide](docs/workflow/README.md) covers all 19 rule files and
16 generated graphs. Each stage guide embeds its graph beside preview and run
commands. Graph arrows show computational dependencies; repeated sample jobs
are collapsed into one box per rule. Some reviewed inputs and standalone
analyses are documented separately, so the graphs do not represent every handoff.

## Data and interpretation

The [sample metadata](figures/cohort/sample_metadata.tsv) and
[detection tables](figures/cohort/) are starting points for exploring the results.
The [input manifest](manifest.csv), [read inventory](data/metadata/raw_files.tsv)
and [source-study records](data/sources/) document the original inputs.
Publication metadata retain updated host names alongside the names used in the
source studies. [Depth conventions](docs/metadata_depth_conventions.md) distinguish
recorded depths from assigned surface values.

All 205 libraries contribute to descriptive summaries. Tests accounting for
sequencing batch use a smaller subset: 136 libraries sequenced on one physical
flowcell, including 122 *Physalia*. Sampling, tissue, host reference availability,
geography and source study overlap. A detection establishes sequence evidence;
it does not by itself establish infection, symbiosis or host specificity.

## Reproducing an analysis

Start with the relevant stage guide above and its
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
