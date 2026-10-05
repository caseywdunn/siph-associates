# siph_associates

Cross-study mining of non-host sequence — **bacteria, viruses, eukaryotic
parasites, and prey** — from whole-genome shotgun data of siphonophores
(Cnidaria: Hydrozoa). This is the public analysis repository for the paper; it
pools three WGS datasets and analyzes them under one symmetric pipeline in which
each study is a first-class covariate, not a hierarchy.

| Study | Data | Role |
|-------|------|------|
| Church et al. 2025, *Curr. Biol.* ([10.1016/j.cub.2025.05.066](https://doi.org/10.1016/j.cub.2025.05.066)) | 151 *Physalia* libraries (5 spp.) + chromosome-scale *P. physalis* genome | Physalia sampling; a host reference |
| Ahuja et al. 2024, *GBE* ([10.1093/gbe/evae048](https://doi.org/10.1093/gbe/evae048)) | 32-species genome skim across the siphonophore phylogeny | phylogenetic breadth |
| Ahuja et al. 2026, *PLoS One* ([10.1371/journal.pone.0351247](https://doi.org/10.1371/journal.pone.0351247)) | *Nanomia* population set + chromosome-scale *N. septata* genome | Nanomia sampling; a host reference |

## Publication preparation (2026-10-03)

Catalog v2, cohort mapping, corrected downstream statistics and eukaryote
reporting have completed validation. Publication figures are generated in
[`figures/`](figures/README.md), with source tables and checksums. The manuscript
and separate supplement are in `../manuscript_siph_associates`.

All 205 libraries remain in descriptive analyses. Only flowcell-based inference
uses the single-physical-flowcell subset (136 libraries overall; 122 Physalia;
12 Nanomia). The genus-wide Nanomia count is 22, including one library from
Ahuja2024 and 21 from Ahuja2026. See the corrected decision records in `docs/`
and the latest RUNLOG entry. Earlier implementation-status sections below are
historical and should not override the accepted versioned outputs.

Collection metadata were enriched on 2026-10-04 from supplied specimen sheets
and 185 public BioSamples. Coordinates now cover all 205 specimens and 169
have tissue descriptions. Depths or intervals cover 202 specimens after
scoring Physalia without an explicit depth record as 0 m, as instructed by
the user. The 150 assigned zeros are distinguished from the 52 source-record
depths by `collection_depth_basis`; see the
[depth convention](docs/metadata_depth_conventions.md).
Existing values and all biological results are preserved. See the
[metadata review](data/metadata/enrichment/2026-10-04/README.md) for source
discrepancies, coverage and reproduction. Manifest rebuilds apply the accepted
`data/metadata/sample_metadata_updates.tsv` overlay.

## Direct trematode comparison and revised figures

A targeted comparison of the eight trematode-positive specimens against
Waki et al.'s deposited markers and their external comparison references
recovers Dinurus-like ITS2-containing and COI sequences and
Prodistomum-like 28S sequences. The comparison retains competing hits and
marker-linkage limits; it does not change the original cohort presence grades.
See [methods, results and reproduction](docs/waki_sequence_comparison.md) and
[reference retrieval](data/sources/waki2026/README.md).

The main-figure set now opens with a study overview and host-photo placeholders,
followed by broad bacterial/archaeal distributions, one focused bacterial
phylogeny, and eukaryotic evidence with host-species labels. Detailed
Mycoplasmatales incidence and gene functions remain available for the
manuscript supplement. See [figure builds and source data](figures/README.md).

## Repository layout

```
README.md              this file
PLAN.md                the analysis plan — symmetric pipeline (S0–S9) + locked parameters
EXECUTION_PLAN.md      goal-ready phase gates, acceptance tests, and SLURM recovery model
WORKFLOW.md            workflow targets, validation, submission, and recovery commands
manifest.csv           unified deduplicated manifest (205 libraries; the S0 input)
data/metadata/         paired FASTQs, provenance, resources, checksums, validation
scripts/
  build_manifest.py    reproducibly builds manifest.csv from data/sources/
workflow/              modular Snakemake rules and atomic rule scripts
config/                production configuration and 205-library sample table
profiles/slurm/         persistent controller/worker submission profile
envs/                   pinned portable software definitions
data/sources/          source supplements from the three studies (provenance)
data/results/          analysis outputs (not tracked)
figures/               publication PDFs, previews, source tables and provenance
```

## The manifest

[`manifest.csv`](manifest.csv) is the single source of truth: **one row per
unique sequencing library** (the analysis unit), with `specimen_id` identifying
the physical animal. It carries taxonomy, collection data, exact depth, explicit
FASTQ lists, true Illumina instrument/run/flowcell/lane batches, and inclusion
state. There are **205 libraries**:
Church 2025 (151 primary rows), Ahuja 2024 (33), and Ahuja 2026 (21). The
study-to-library relationship is normalized separately in
[`data/metadata/library_provenance.tsv`](data/metadata/library_provenance.tsv):
Church YPM-IZ-104465/Ahuja NA22 is one library with two provenance records;
CWD16 and NA19 are distinct Ahuja 2024 libraries both reused by Ahuja 2026.

Explicit selected mate pairs and sequencing batches parsed from FASTQ headers are in
[`data/metadata/raw_files.tsv`](data/metadata/raw_files.tsv). Phase-0 resources,
host-reference audits, validation, and the checksum freeze are in the same
directory. Twenty-one of the 151 published Church libraries were submitted to
PRJNA1092115 but are **not yet released by NCBI**, so
their `sra_run` is left blank (reads held locally; release pending). The only
field not fully resolved is per-specimen *Nanomia* CO1 species for the 7
Rhode-Island / Hawai'i congeners (marked `Nanomia sp.`; needs Ahuja's CO1 calls).

Rebuild it with:

```bash
python3 scripts/build_manifest.py
python3 scripts/inventory_phase0_resources.py
python3 scripts/summarize_phase0.py
python3 scripts/validate_manifest.py
python3 scripts/freeze_manifest.py
python3 scripts/validate_manifest.py --freeze data/metadata/manifest.freeze.sha256
```

## Status

Phases 0--2 and the Phase-3 strategy benchmark are complete. The 205-library manifest and normalized input
metadata are validated and checksum-frozen, and all 205 libraries passed the
universal fastp, Kraken2/Bracken, sylph, and phyloFlash screen. The accepted
Phase-2 aggregation contains 2,282,231 nomination-only rows; these are not
presence claims. Phase 3 is active: the accepted benchmark selects reference
depletion for reference-bearing hosts and deterministic 25 M-pair trimmed-read
inputs for reference-free hosts. All 205 cohort inputs are prepared; 195
libraries meet the 175,000-pair assembly floor and are frozen in
`config/phase3_assembly.tsv`. The per-library assembly stage (`phase3_assembly`)
is implemented and fixture-tested; see `WORKFLOW.md`.
The manuscript (LaTeX) lives in a separate repository and cites this one.

## Citation / license

TBD (add Zenodo DOI on release).
