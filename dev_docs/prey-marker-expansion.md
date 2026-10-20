# Fish host-depleted expansion and branch separation

## Final reporting assessment, 2026-10-19

Array 13512346 completed all 44 tasks with exit code 0. Final reporting
reclassified the 161 host-depleted libraries with accession-linked reference
metadata and combined all 205 libraries without new read searches. All report
hashes, cohort membership, unique read-pair identities and category totals passed.
There are no identical paired sequences among the 32 strong matches. Molecular
independence and sequencing contamination remain unresolved; no species-level
prey assignments or counts of prey individuals are asserted.

The [final report](../docs/fish_12s_results.md) supersedes the running status
and preliminary genus-based classification below. Generated results and source
hashes are in `figures/fish_12s/`. The completed controls and integration fixture
below remain applicable because sequence-search settings are unchanged.

## Full-read fish expansion, 2026-10-18

Casey approved screening all 44 remaining libraries at up to 200 million raw
pairs per library, superseding the proposed diet-selected subset. This adds
7,939,126,078 pairs before QC. Crustacean work remains paused.

`config/fish_12s_full_reads.tsv` freezes the raw paths and lane order from the
original sequence-analysis inventory. `workflow/fish_12s_full_reads.smk` reuses
the original `trim_reads.py` and its fastp/BBTools environment, writing only to
the new fish run directory. It shares the unchanged filter commands and BLAST
settings with the host-depleted workflow. The regenerated inputs must match
previous input hashes, settings and post-QC counts before screening.

The previous control checker inspected recruitment but not BLAST evidence.
It now uses the same competitive classifier as the reports and verifies four
fish positives and exclusion of the frog. Missing/empty alignments and a frog
misclassified as fish fail regression tests. Production launchers re-evaluate
the controls rather than trusting the old marker.

RefSeq classification now uses accession-linked GenBank lineage annotations
instead of a MitoFish-derived genus list. Reference job 13512182 exposed
CONTIG-only GenBank records and failed without producing a metadata table;
job 13512217 completed after adding explicit accession-version/length checks
for those records. All records with sequence content were checked against the
frozen FASTA. Metadata, exception accessions and source snapshots are under
`data/results/fish_12s/reference_metadata/2026-10-18/`.

Validation job 13512256 exercises the production path on two lanes with both a
forced exact cap and an uncapped input. It completed successfully in 1:50:
60 pairs entered the capped screen (33 candidates), and all 90 entered the
uncapped screen (50 candidates). Both retained all represented fish controls
and rejected the frog at competitive classification. Ten regression tests,
Ruff, shell syntax and both generated graph freshness checks passed.
The full production dry run schedules
44 preparations and 44 instances of each screening/review stage (309 steps
including the aggregate target). The host-depleted dry run reports nothing to
do. The original analysis Snakefile and sequence rules are unchanged.

Production array **13512346** was submitted with
`sbatch --array=0-43%4 scripts/fish_12s/batch_full_reads.sh`. Tasks 0–3 were
confirmed running and preparing reads; tasks 4–43 wait on the concurrency limit.
The first four samples are Ahuja2024__CWD146 (Erenna), Ahuja2024__CWD131 and
Ahuja2024__CWD20 (Rhizophysa), and Ahuja2024__CWD11 (Stephanomia). This ordering
only prioritizes execution; all 44 libraries remain selected. Source snapshots,
input identities and completion status are recorded separately for each task.
The cohort is **running, not yet complete**. Keep computational source and
configuration stable while the array runs. Logs are
`logs/fish_12s/full-reads-13512346_{0..43}.log`.

The [execution guide](../docs/prey_marker_screen.md#full-read-fish-expansion)
contains the canonical launch commands and output layout. The 161-library
tables below retain their original preliminary classifications; accession-based
reclassification is a reporting step and does not require new sequence searches.

## Completed fish cohort

The [cohort manifest](../config/prey_marker_host_depleted.tsv) contains all 151
Physalia and ten Nanomia septata libraries: 416,741,747 retained read pairs.
There is no additional subsampling. The other 44 libraries are not yet screened.
The [execution guide](../docs/prey_marker_screen.md#all-host-depleted-libraries)
provides launch, control and dry-run commands.

Fish job 13501886 completed all 625 new-workflow steps in 1:59:40. It reused five
completed pilots and screened 156 additional libraries. Control job 13501764
passed the same cohort implementation on nine synthetic pairs. The
[cohort review](fish-12s/positive_pilot.md) preserves the 161-library summary and
all 1,021 candidate pairs; 32 pairs in eight Physalia libraries have long, close
fish matches under the preliminary review criteria.

The historical submissions and frozen executed source remain in
`data/results/prey_marker_submissions/2026-10-18-verified-expansion/`.
`submission.json` records job identities and the submitted archive hash.
Per-job commands, software, input identities and outputs remain under
`data/results/fish_12s/executions/`. The original fish submission waited for the
then-shared COI verification job 13501446. That historical scheduling decision
is preserved in the submission record, but COI verification is no longer a
prerequisite for running fish analyses. Fish controls remain required.

## Crustacean investigation paused, 2026-10-18

At Casey's request, crustacean source, reference audits, controls and review
records moved to branch `crustacean-coi-screen`. That branch preserves
`dev_docs/crustacean-coi/verification.md` with the completed verification results
and resumption considerations. Cohort job 13501881 was cancelled after 8:58:43;
its last logged progress was 93/1240 steps. Partial outputs are not a completed
cohort result. All completed and partial data remain under
`data/results/crustacean_coi/`, outside Git; nothing is resubmitted.

The active `fish-12s-screen` branch contains the fish code and reports. Its
launcher no longer reads a COI gate or archives COI code/configuration. The
manifest generator now writes only fish manifests, and the graph registry
contains only the fish and established manuscript workflows.

The existing cohort manifest keeps its historical shared-panel selection label
and byte-identical content. The fish search scripts, settings and workflow rules
are unchanged. When Git checkout restored identical files, their pre-separation
modification times were retained to avoid invalidation from checkout alone.
Sequence outputs, metadata identifiers and the accepted manuscript workflow are
unchanged. The separate crustacean branch preserves the original combined plan
and all evidence needed to revisit it.

Separation checks passed: all 75 pre-separation changed/untracked files are
preserved on the crustacean branch (apart from explicit pause notes and one
trailing-blank-line cleanup). Fish workflow, search/report scripts, cohort
configuration, frozen manifest and result tables remain byte-identical. Both
manifest-generation modes retain the expected members, counts and resolved
paths, and neither writes a crustacean manifest. Production and synthetic fish
launcher checks pass without COI files. Ruff, shell syntax, four fish-reference
tests, the fish rule-graph freshness check, and local documentation links pass;
16 crustacean checks passed before archiving. The production fish dry run reports
"Nothing to be done (all requested files are present and up to date)." No sequence
analyses ran during separation. The dry-run record is retained under
`data/results/branch_separation/2026-10-18/`.
