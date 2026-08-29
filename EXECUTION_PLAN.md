# Publication-analysis execution plan

Last updated: 2026-08-29.

This document turns the scientific design in `PLAN.md` into a sequence that can
be pursued as one persistent goal. The exploratory repository
`../20260513_siph_symbionts` is evidence for decisions and a source of tested
implementation patterns; it is frozen as a pilot. Publication outputs are
recomputed in this repository from the unified manifest and a versioned
configuration.

## Goal and completion criteria

**Objective:** build, run, validate, and document a restartable, study-agnostic
pipeline that applies the same read preprocessing and detection screens to all
unique siphonophore libraries, constructs pooled associate catalogs, validates
candidate presences by mapping breadth or assembled sequence, and produces the
tables, figures, provenance, and Methods-ready records needed for publication.

The goal is complete only when:

1. the analysis manifest contains one row per unique raw library, has no
   unresolved paths or accidental duplicate specimens, and passes automated
   validation;
2. every production parameter, database, reference, software environment, and
   decision rule is recorded before the corresponding full-cohort run;
3. all expected per-library and cohort outputs exist, pass integrity checks,
   and can be recreated by the workflow without relying on the exploratory
   output tree;
4. bacterial/archaeal, viral, parasite, and prey candidates have explicit
   evidence grades; k-mer calls alone are never reported as presence;
5. primary results, sensitivity analyses, figures, and Methods/provenance are
   generated from machine-readable publication-repository outputs; and
6. a clean dry-run finds no missing work and the final validation report has no
   unexplained failures or exclusions.

## Decisions to lock before cohort computation

These are gates, not details to decide opportunistically after seeing the final
results.

### D0.1 Analysis units and manifest

Phase 0 resolved the identities to **205 unique raw libraries**. Church
`YPM-IZ-104465` and Ahuja 2024 `NA22` point to the same raw library and are one
analytical row with two provenance records. NA19 and CWD16 are distinct Ahuja
2024 libraries that are both explicitly reused in Ahuja 2026 Table S8; each is
one analytical row with two study-provenance records. NA19 was additionally
delivered in both local data trees.

Define three identifiers explicitly:

- `library_id`: unique sequencing library/read set and workflow wildcard;
- `specimen_id`: physical animal and statistical sampling unit;
- `provenance_record`: study-specific occurrence of that library/specimen.

The manifest now carries `sequencing_batch`, explicit raw R1/R2 lists, read
length, `include_primary`, and `exclusion_reason`; normalized raw-file and
study-provenance tables retain the one-to-many detail. The stale paths for
Church libraries YPM-IZ-111839 and
YPM-IZ-111840 were corrected on 2026-08-29 to their extant locations under
`illumina/physalia_novaseq9_2023/`; both mates and the stored read-pair counts
were verified. The production tables are checksum-frozen and the automated
validator checks unique IDs, paired files, readability, pair counts, categories,
reference routing, resources, and host-audit evidence.

### D0.2 Scope and strength of claims

The primary product is a breadth-validated catalog of detectable associates,
not relative microbial abundance: host DNA fraction and sequencing depth vary
too much for abundance comparisons. `read_pairs`, study, and sequencing batch
are detection-sensitivity covariates. Cross-host patterns are descriptive where
host species have one specimen; inferential host-species/geographic models are
restricted to replicated strata (principally *Physalia* and *Nanomia*) and use
study/batch-aware or restricted permutations.

Use evidence grades consistently:

- **nominated:** Kraken2/Bracken, sylph, phyloFlash, or marker-search evidence;
- **validated:** competitive mapping clears a predeclared breadth/read support
  rule, or an assembled sequence has concordant taxonomic evidence;
- **high confidence:** independent mapping plus assembly/marker evidence;
- **trace/ambiguous:** below the validation rule or consistent with conserved
  sequence/repeat recruitment;
- **probable contaminant:** breadth and cross-library/batch pattern support
  contamination; organism identity alone is insufficient.

Before full mapping, lock a bacterial presence rule from pilot positives,
negative/decoy references, and mapping nulls. A reasonable starting rule to
evaluate is >=10% genome breadth with reads filtered at >=95% aligned identity,
MAPQ >=30, proper pairs, and a minimum pair count; report breadth continuously
and run threshold sensitivity analyses. Do not confuse CoverM's per-read
`--min-covered-percent 75` with 75% genome breadth.

### D0.3 Reference-free assembly strategy

Do not label Kraken2 `--classified-out` as general “non-host” data. It is a
taxonomically nominated subset and will omit novel bacteria, viruses, and
eukaryotes absent from the database. Before cohort assembly, benchmark on a
small, predeclared panel spanning high/low signal and both reference states:

1. whole-genome reference depletion where a conspecific reference exists;
2. Kraken taxon-nominated read assembly; and
3. fixed-effort assembly of unfiltered trimmed reads for reference-free hosts,
   followed by contig-level host/non-host classification.

Compare non-host marker recovery, viral recovery, MAG recovery, host carryover,
runtime, memory, and disk. Lock a feasible method for the 44 reference-free
libraries. If fixed-effort whole-read assembly is selected, use a deterministic
paired subsample and state its detection limit. The assembly branch is a
reference-building/validation branch; universal read screens and final mapping
remain the cross-library comparison, so unequal assembly yield cannot become an
abundance measure.

### D0.4 Domain-specific rules

- **Bacteria/archaea:** pooled nomination; taxon-verified NCBI references plus
  dereplicated MAGs; one competitive catalog mapped to every library.
- **Viruses:** lock geNomad/CheckV quality criteria, vOTU dereplication identity
  and coverage, and whether host linkage uses iPHoP, CRISPR spacers, or both.
  Host-reference EVEs are reported separately from associate viruses.
- **Parasites/prey:** lock PR2 version and non-host rules; use 18S/28S and COI or
  mitogenome evidence. Parasite versus prey is a biological interpretation with
  an explicit evidence field, not a classifier taxon bin. Targeted mapping is
  only used with a sufficiently close reference and a marker-appropriate
  breadth rule.
- **Contamination:** derive batch from raw file/run metadata; retain common
  contaminant taxa in QC outputs and test their breadth/batch distribution.

## Implementation and execution phases

Each phase ends with a committed checkpoint and an updated `RUNLOG.md` giving
the workflow target, controller job ID, worker-job query, output location,
validation result, and exact next command. Do not begin a full-cohort phase
until its gate is accepted.

### Phase 0 — reconcile and freeze inputs

**Status: complete (2026-08-29).** The checksum-aware validator passes all
acceptance criteria for 205 analytical libraries, 313 paired FASTQ records, 208
provenance records, 10 locked resources, and both host-reference audits. No
cohort analysis was started. See `RUNLOG.md` and
`data/metadata/manifest_validation.txt`.

Deliverables:

- corrected `manifest.csv` plus normalized provenance and raw-file tables;
- `scripts/validate_manifest.py` and a checked validation report;
- checksummed host references, databases, and source metadata;
- resolved *Nanomia* species calls where source data permit, otherwise explicit
  unresolved labels that are not silently used as species-level observations;
- reference-scaffold audit summaries imported as provenance, with the already
  completed *P. physalis* and *N. septata* pilot findings independently checked;
- a storage estimate based on the manifest, not a hard-coded library count.

Acceptance: zero duplicate raw-file sets, zero missing R1/R2 inputs, zero
unexplained manifest TODOs, and the expected number of unique analysis rows.

### Phase 1 — build the reproducible workflow skeleton

Implement a Snakemake workflow in the repository root with modular rules,
`config/config.yaml`, `config/samples.tsv`, pinned environment definitions, and
a SLURM profile. Port small tested functions from the pilot only after making
them sample-table-driven and removing study-specific paths.

Required engineering properties:

- atomic outputs followed by format/content validation;
- per-rule logs, benchmarks, software versions, parameter snapshots, and input
  checksums;
- `--rerun-incomplete`, `--keep-going`, a generous filesystem latency wait, and
  explicit retries only for known transient failures;
- persistent, regenerable trimmed reads on scratch; publication products on
  `/work`; no dependence on files under the exploratory repository;
- a `rule all` plus stage targets and validation sentinels;
- dry-run, DAG, lint, and small fixture/smoke tests before real data;
- no silent empty-output success and no output shared by concurrent array jobs.

Scratch files are normal workflow outputs, not the scientific record. Missing
scratch inputs are regenerated only when pending downstream work requires them;
completed persistent products must not be invalidated merely by scratch purge.
Record forced reruns when code or parameters change under mtime-based triggering.

Acceptance: lint and dry-run pass; a representative library from each study and
host-routing class completes trim, all three screens, host handling, and output
validation end to end.

### Phase 2 — universal preprocessing and screening

For all included unique libraries:

1. deterministically cap only libraries above 400 M raw read pairs, preserving
   mate synchronization and recording seed and retained fraction;
2. run fastp with the locked settings;
3. run Kraken2/Bracken, sylph against GTDB r220 plus OceanDNA, and phyloFlash
   against the locked SILVA build on the same trimmed reads; and
4. aggregate library QC and candidate nominations without making presence
   claims.

Run a 6–10 library resource pilot first, then set rule-specific SLURM resources
from observed maxima with margin. Validate paired-read counts, fastp reports,
classifier completion, expected sample columns, and aggregation completeness.

Acceptance: one validated output set per included library for every screen, with
all failures either rerun successfully or explicitly excluded with a reason.

### Phase 3 — lock and run the assembly branch

Run the D0.3 benchmark, commit its decision table, then execute host handling and
assembly. Reference-bearing hosts map to the complete audited nuclear plus
mitochondrial reference and retain both-unmapped pairs. Reference-free hosts use
the benchmark-selected approach. Set the assembly inclusion rule from observed
read yields before selecting libraries; never hand-pick them from interesting
taxonomic results.

Run MEGAHIT, contig QC, mapping back where binning requires coverage, MetaBAT2,
CheckM2, GTDB-Tk, geNomad, CheckV, barrnap/marker extraction, and cohort
dereplication. Preserve per-library provenance for every MAG, vOTU, and marker.

Acceptance: every eligible library has a validated assembly or a documented
objective exclusion; all retained MAGs/vOTUs/markers meet locked QC criteria and
carry source-library and method provenance.

### Phase 4 — build frozen pooled catalogs

Construct catalogs without study-specific tiers:

- nominate across all Phase-2 screens;
- add taxon-verified fetched genomes and Phase-3 MAGs;
- dereplicate and screen for contaminated/mislabeled references;
- build viral vOTUs and eukaryotic marker/reference sets under the locked rules;
- assign stable catalog IDs and write a catalog manifest containing accession,
  version, taxon, sequence checksum, evidence, source, and exclusions.

Freeze each catalog before cohort mapping. Later additions create a new catalog
version and require a complete remap, never an incremental subset.

Acceptance: catalog validation passes, every sequence is attributable, and
positive/negative/decoy mapping tests support the chosen presence rules.

### Phase 5 — uniform validation and quantification

Map every library against each applicable frozen catalog using identical
filters. Produce long-form tables of breadth, read pairs, aligned identity,
depth, catalog version, and evidence grade. Quantify viruses by competitive
mapping to the dereplicated vOTUs; quantify eukaryotic targets only where
reference divergence and marker structure make breadth interpretable.

Acceptance: complete library-by-reference matrices, no duplicated sample/catalog
keys, mapping totals reconcile with per-library QC, and pilot headline findings
are treated as hypotheses that may be confirmed or overturned.

### Phase 6 — statistical analysis and robustness

Predefine primary and sensitivity analyses before inspecting final plots:

- catalog incidence and evidence-grade summaries by host and study;
- replicated-stratum models for *Physalia* and *Nanomia*, adjusting for depth,
  geography, study, and sequencing batch as identifiable;
- host-phylogeny summaries across the broad skim with explicit acknowledgment
  of singleton species and study confounding;
- reference-based versus reference-free host-handling concordance for the two
  reference species;
- presence-threshold sensitivity and capped-versus-uncapped checks on a
  stratified subset;
- contamination/batch tests and leave-one-study-out catalog sensitivity; and
- cross-domain co-occurrence described as association, with phage host links
  supported by sequence-based linkage rather than co-occurrence alone.

Acceptance: analysis scripts regenerate all reported tables/figures from frozen
long-form results; claims are stable or their sensitivity is reported.

### Phase 7 — publication freeze and reproducibility audit

Generate Methods-ready parameter/version tables, sample inclusion flow,
database/reference citations, QC summaries, figures, and machine-readable
source-data tables. Run the workflow from a clean metadata state in dry-run
mode, audit expected versus observed artifacts, tag the catalog and analysis
release, and archive checksums plus the final config. Large reads and databases
remain external and are documented by stable accessions/paths and checksums.

Acceptance: `snakemake -n` reports nothing to do for publication targets; the
validation report is clean; every manuscript number has a producing script and
source table.

## SLURM and interrupted-session operating model

Never run the Snakemake scheduler itself only inside the interactive tmux
allocation. Submit a lightweight **controller job** with `sbatch`; it launches
worker jobs through the SLURM profile. Both controller and workers then survive
loss of the interactive session. Use stage targets so a controller fits the
available partition wall time. If it times out, completed outputs remain and the
same target is safely resubmitted with `--rerun-incomplete`.

Operational rules:

1. Before submission, commit code/config and update `RUNLOG.md` with the commit,
   target, planned resources, and command.
2. After `sbatch`, immediately record the controller job ID. Worker IDs are
   recoverable by job name and Snakemake logs.
3. Use dependencies for explicit one-off jobs; let Snakemake express ordinary
   DAG dependencies.
4. On reconnect, inspect `RUNLOG.md`, `squeue`, `sacct`, workflow logs, and the
   target validation report before resubmitting anything.
5. Never infer success only from SLURM state: require validated output sentinels.
6. Keep no irreplaceable state in tmux, local shell variables, or dotfiles that
   are not committed/documented.

## Recommended goal-command objective

Use the objective at the top of this document, with these phase gates. The first
goal turn should execute **Phase 0 only**, update the plan with resolved counts
and decisions, and stop at the D0.3 benchmark design if a scientific choice is
still needed. Subsequent turns should resume from `RUNLOG.md` and the workflow
validation state, not from conversational memory.
