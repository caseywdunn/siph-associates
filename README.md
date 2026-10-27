# Organisms associated with siphonophores

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23441461.svg)](https://doi.org/10.5281/zenodo.23441461)

This repository contains the analyses accompanying:

**Organisms associated with siphonophores revealed by short-read whole-genome sequencing**

*Manuscript citation and public link to be added.*

We examine DNA from bacteria, viruses, prey and parasites captured
alongside siphonophore DNA during genome sequencing. The study screens 213 specimens:
151 *Physalia*, 22 *Nanomia* and 40 other siphonophores. Biological summaries
use 209 specimens after four library exclusions and set aside the specified
contaminant taxa. All 151 *Physalia* libraries remain in those summaries. Start with the tables
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
table identifiers refer to the manuscript and S1 Appendix. Figures A–H and
Tables A–B are internal items in that appendix.
Each guide explains the methods and links to their implementation; the data links
point to small tables retained in this repository.

| Results section in the paper | Figures and tables | Analyses and files |
| --- | --- | --- |
| Human and laboratory-associated sequences | Fig H in S1 Appendix | [Human-read methods](docs/human_read_screen.md), [excluded libraries](config/excluded_libraries.tsv), [contaminant taxa](config/contaminant_taxa.tsv), and [human/C. acnes distributions](figures/contamination/human_cutibacterium.pdf) |
| Genome sequencing reveals a broad assemblage of associated organisms | Fig 1; Figs A–B in S1 Appendix; Table 1 | [Workflow overview](docs/workflow/README.md), [sample metadata](figures/cohort/sample_metadata.tsv), [sampling](figures/cohort/host_sampling.tsv), and [evidence summary](figures/overview/domain_summary.tsv) |
| Bacterial associates include diverse proteobacteria and recurrent Metamycoplasmataceae | Fig 2; Fig E in S1 Appendix; Table 1 | [Genome recovery](docs/workflow/phase3_assembly/README.md), [genome selection](docs/workflow/phase3_catalog/README.md), [detection criteria](docs/workflow/phase6_analysis/README.md), and [host distributions](figures/overview/microbial_host_distribution.tsv) |
| *Vibrio* diversity extends across *Physalia* and *Nanomia* | Text results | [Read mapping](docs/workflow/phase5_mapping/README.md), [assembly support](docs/workflow/phase6_analysis/README.md), and [bacterial detections](figures/cohort/bacterial_detections.tsv) |
| Siphonophore Metamycoplasmataceae have close relatives in other marine animals | Fig 3; Figs C and E in S1 Appendix | [Genome phylogeny](docs/workflow/mycoplasmatales_phylogeny/README.md), [cnidarian marker comparisons](docs/workflow/cnidarian_16s/README.md), [Cassiopea comparisons](docs/workflow/cassiopea_16s/README.md), [tree data](figures/mycoplasmatales/source_phylogeny_focus.tsv), and [host incidence](figures/cohort/mycoplasmatales_incidence.tsv) |
| Bacterial communities differ among geographic regions and hosts | Figs D and F; Table A in S1 Appendix | [Regional tests](docs/workflow/phase6_analysis/README.md), [matched host and regional comparisons](docs/coauthor_metadata_analysis.md), [incidence counts](figures/coauthor_metadata/targeted_incidence.tsv) |
| Uneven sampling limits comparisons among tissues and depths | Figs C–D and G in S1 Appendix | [Tissue and depth methods](docs/coauthor_metadata_analysis.md), [tissue incidence](figures/coauthor_metadata/tissue_incidence.tsv), [depth design](figures/coauthor_metadata/host_depth_design.tsv), and [fish detections](figures/fish_12s/sample_summary.tsv) |
| Viral sequences include bacteriophages and endogenous candidates | Table 1 | [Viral recovery and selection](docs/workflow/phase3_catalog/README.md), [viral catalogs](docs/workflow/phase4_catalog/README.md), and [viral detections](figures/cohort/viral_detections.tsv) |
| Eukaryotic associates include prey and parasite lineages | Fig 4; Table 1 | [Ribosomal-marker analysis](docs/workflow/eukaryote_gate/README.md), [eukaryote evidence](figures/cohort/eukaryote_evidence.tsv), [fish 12S results](docs/fish_12s_results.md), [trematode comparison](docs/waki_sequence_comparison.md), and [Perkinsela bin review](docs/mag_domain_audit.md) |
| Assembly recovered bacteria absent from reference databases | Text and supplementary results | [Screening methods](docs/workflow/screens/README.md) and [screen comparisons](figures/screening/) |
| Host DNA removal changed some detections | Text and supplementary results | [Host-reference comparisons](docs/coauthor_metadata_analysis.md#host-reference-context), [current removal statistics](figures/coauthor_metadata/host_depletion_summary.tsv), and [historical fixed-catalog comparison](figures/coauthor_metadata/historical_host_concordance.tsv) |

Table B in S1 Appendix inventories the accompanying data.
The [reviewed catalog](figures/manuscript_revision/catalog/mag_catalog.tsv) provides genome quality.
The [supporting-data tables](tables/README.md) provide phylogenetic and marker
comparisons, statistical results, and trematode-marker evidence
with original paths and checksums.
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
| Assess reads matching the human genome | Screen bounded raw-read prefixes with stringent paired alignments, simulated controls and candidate comparisons against broader references | [Methods and interpretation](docs/human_read_screen.md), [validation](dev_docs/human-screen/2026-10-22-validation.md) |
| Compare bacterial relationships | Infer a genome tree and search earlier marker surveys | [Phylogeny](docs/workflow/mycoplasmatales_phylogeny/README.md), [tree graph](docs/workflow/graphs/phylogeny.svg), [16S comparisons](docs/workflow/cnidarian_16s/README.md) |
| Compare distributions and assess sensitivity | Describe tissue and depth patterns, compare regions within hosts and hosts within regions, and assess sequencing-batch effects | [Statistical methods](docs/workflow/phase6_analysis/README.md), [current reporting rules](workflow/manuscript_revision.smk), [collection context](docs/collection_context_analysis.md), [matched comparisons](docs/coauthor_metadata_analysis.md) |
| Produce figures and supporting tables | Summarize accepted results for the paper | [Current reporting guide](docs/expanded_reporting.md), [graph](docs/workflow/graphs/manuscript_revision.svg), [figure and data index](figures/README.md) |

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

Full evidence tables retain all 213 screened libraries and flag the four
excluded libraries and contaminant taxa. Biological summaries use the 209
retained specimens. Tests with permutations within sequencing flowcells use
libraries sequenced on one physical flowcell; the generated
[flowcell table](figures/cohort/library_flowcells.tsv) records eligibility. The
Physalia geographic tests retain their 122 eligible libraries. Sampling, tissue, host reference availability,
geography and source study overlap. A detection establishes sequence evidence;
it does not by itself establish infection, symbiosis or host specificity.

## Reproducing an analysis

Start with [current cohort and reporting](docs/expanded_reporting.md), which runs
[`workflow/manuscript_revision.smk`](workflow/manuscript_revision.smk) through its
[preview and launch commands](docs/expanded_reporting.md#preview-and-launch),
then the relevant stage guide and its
[environment and input instructions](docs/workflow/README.md#environment-and-inputs).
The workflow uses Snakemake to connect analysis steps and their files. Each guide
provides a separate **dry run** (a preview that executes no analysis jobs), a
launch command and a generated rule graph. Instructions for Yale's computing
cluster and job recovery are in the [operation guide](WORKFLOW.md).
The default Snakemake target checks input readiness; it does not run the whole paper.

This repository retains code, settings, small source files, figures and supporting
tables. Large sequencing files, reference databases and most sequence-analysis
outputs are stored separately; accessions and configured locations are described
in the input instructions. Public run accessions are recorded for 211 of the 213 libraries; OA-190531-1
and OA-Bluebottle still lack verified accessions. The repository includes selected
trematode-marker sequences, but a complete public deposit of reconstructed genomes
and markers remains pending. Some steps also require reviewed intermediate inputs.
A clone alone is therefore insufficient for a complete rerun. The stage guides identify
these requirements and where outputs are written.

The [NCBI submission guide](docs/ncbi_submission.md) provides a progress
checklist and instructions for BioProject/BioSample registration, bacterial
genomes, viral and marker sequences, validation, accession tracking and release.
It identifies the remaining metadata and submission-route decisions for these
reused sequencing data. The [sequence preparation guide](docs/deposit_preparation.md)
and [local validation guide](docs/deposit_validation.md) describe the existing
cluster commands and draft outputs. The provisional
[organism-name crosswalk](config/mag_organism_names.tsv) records proposed NCBI
names and unresolved entries. Deposit preparation preserves the mapping
catalogs; catalog read mapping is not repeated. Submission and publication of
the sequence archive remain pending.

## Citation and development

The v0.5.0 code release is archived on Zenodo at
[doi:10.5281/zenodo.23441461](https://doi.org/10.5281/zenodo.23441461).
That archived release precedes the current 213-library reporting and does
not contain the reconstructed sequence deposit. Use the files linked above for
the current paper; the archive is a versioned code snapshot.
The manuscript citation and public supplementary-material link will be added when
available. [Workflow documentation checks](dev_docs/workflow-documentation.md)
and [execution records](RUNLOG.md) are maintained separately from the methods guides.
