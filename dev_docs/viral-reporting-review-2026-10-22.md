# Coauthor reporting review, 22 October 2026

The bounded workflow `workflow/coauthor_reporting.smk` applied the reviewed
Negarnaviricota interpretation and the publication metadata's identification of
NA37 as *Nanomia bijuga*. It read completed v3 evidence matrices and statistical
outputs. No FASTQs, alignments or assemblies were rebuilt.

The dry run contained eleven jobs: metadata, viral classification, three
flowcell-test jobs, three incidence jobs, validation, publication exports and
the aggregate target. The final dry run reported nothing to do. Evidence
validation passed with zero errors; two focused classification tests and Ruff
checks passed. The no-change check verified all mapping measurements, evidence
grades, target IDs and sequence checksums. Existing Physalia regional outputs
were reused only after verifying identical Physalia metadata and bacterial
evidence. Host-stratified flowcell tests and incidence were recomputed.

The derived catalog has 404 associate and 160 endogenous-candidate groups.
Seven accepted VOTU00242 observations, including three with assembly support,
are now reported as possible endogenous elements. VOTU00144 and VOTU00165 have
no accepted observations. Associate totals are 624 detections, 391 detected
groups, 122 specimens and 439 assembly-supported detections. Excluding the
background-flagged formalin specimen leaves 538 viral detections, including
353 with assembly support. All bacterial and eukaryotic detection counts are
unchanged.

There are now 202 specimens with named host identifications and eleven with
unresolved identifications. None of the three bacterial targets in the pooled
Nanomia flowcell comparison is testable after the NA37 correction: one lacks
sufficient eligible detections and two have constant statistics. No specimen
was reprocessed. Fish report host names were also refreshed from the current
publication metadata without rerunning its screen.

Life-stage summary generation and its table/macros were removed because the
labels are inconsistent. Raw metadata values remain available.

The current reporting selector points to
`data/results/coauthor_review/2026-10-22/phase6_analysis/v3` and
`figures/coauthor_review/catalog`. The original 21 October statistics and v3
mapping catalogs are preserved. In the run directory, `prepared.json` records
source hashes and reused-output hashes, `validation.json` records refreshed
table hashes and the validation command, `execution.json` records commands and
source hashes, and `provenance/` contains executed code and configuration
snapshots. `dry_run.log`, `statistics_run.log`, `publication_run.log` and
`final_dry_run.log` retain execution evidence.

Public rebuild instructions and the generated rule graph are in
[`docs/viral_reporting_review.md`](../docs/viral_reporting_review.md).
