# COI origin verification, 2026-10-18 UTC

This follow-up implements the next step in [the candidate audit](validation.md).
Methods, thresholds and reproduction commands are in the
[analysis guide](../../docs/crustacean_coi_pilot.md#verification-with-prokaryotic-competitors).
The investigation began on `fish-12s-screen`; its current source and review
records are now preserved on `crustacean-coi-screen`.

## Paused on a separate branch, 2026-10-18

Casey paused the crustacean investigation. Its source, controls, reference audits,
and review records are preserved on `crustacean-coi-screen`; fish development
continues on `fish-12s-screen` independently.

Verification job 13501446 completed successfully in 1:18:20 and passed its gate.
Of the 92 candidate pairs from 2,152,939 host-depleted pairs in ten *N. septata*
libraries, 39 were bacterial/archaeal oxidases, one supported another eukaryote,
15 were ambiguous, and 37 remained unresolved. None supported a crustacean call.
The two NA35 candidates were bacterial false positives. All 115 divergent
positive controls were retained, but only three supported a crustacean call;
the absence of supported real calls is not evidence of prey absence. The
[completed summary](validation/verification_summary.json) and
[passed gate](validation/verification_gate.json) preserve these results in Git.

Cohort job 13501881 was cancelled at Casey's request after 8:58:43. Its last
logged progress was 93/1240 steps. Scheduler accounting records cancellation;
the launcher's `exit_code.txt` says zero because its EXIT trap did not propagate
the cancellation, so that file does not establish completion. Partial cohort
outputs are not a completed result. The job, execution log, frozen submission
source, completed pilot and verification outputs remain under
`data/results/crustacean_coi/` and
`data/results/prey_marker_submissions/2026-10-18-verified-expansion/`.
These bulk outputs are ignored by Git and remain in this working directory.

Before resuming, review the conservative assignment sensitivity and distinguish
completed outputs from interrupted files using scheduler/execution records and
Snakemake metadata. Preview the cohort workflow on this branch before any new
submission; no resumption is scheduled. The original assembly, catalog, mapping,
and manuscript analyses are unchanged.

The new reference set contains 42,947 UniProt bacterial/archaeal oxidases
(release 2026_03, released 2026-09-02), plus 17,781 annotated RefSeq mitochondrial
proteins. Lineage review assigns 669 mitochondrial proteins to crustaceans and
17,112 to other eukaryotes. It corrects 139 nucleotide-reference records labeled
non-target solely because they fell outside the original seven target groups.
The original panel and pilot outputs are preserved.

All 108,778 target proteins are audited. Exactly identical proteins are searched
once (21,447 unique sequences), then the evidence is transferred to every source
record. Source-specific manual quarantine still applies independently. Automatic
quarantine identifies unreliable assignment labels; it does not establish that
every quarantined record is a bacterial contaminant.

The output root is `data/results/crustacean_coi/verification/2026-10-18/`.
`source/` preserves the UniProt response and release identity, `anchors/` the
lineage correction and reference audit, `curated/` the assignment panels, and
`verification/` the candidate/control alignments. `summary.json`,
`classified_pairs.tsv` and `gate.json` are written only after the searches finish.

| Execution | State when recorded | Scope |
| --- | --- | --- |
| 13501245 | Cancelled after 4 min | Replaced redundant searches of identical target proteins |
| 13501353 | Failed before computation | Stale Snakemake lock; cleared after both jobs were confirmed stopped |
| 13501446 | Running | Audit unique reference proteins, then verify all 92 real pairs and both control sets |

Each execution has a recoverable source snapshot and command/resource records
under `data/results/crustacean_coi/executions/{job}/`. The original real-library
FASTQs are not scanned in this verification. Original recruitment remains broad;
all candidates are preserved, including ambiguous and divergent matches.

Read `gate.json` and scheduler state before interpreting this run as complete.
The gate requires both known NA35 bacterial pairs to be recognized, all four
flagged barcodes quarantined, all 92 pilot pairs accounted for, no supported
crustacean calls among negatives, retention of the 115 divergent positive controls,
and at least one supported original positive control. Retaining a divergent
candidate is distinct from making a confident taxonomic assignment.

The conditional cohort jobs and their frozen source are recorded in
[the expansion notes](../prey-marker-expansion.md). No results from this work have
been incorporated into the manuscript.
