# Fish host-depleted expansion and branch separation

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
