# Host-depleted marker expansion, 2026-10-18 UTC

[The methods and commands](../docs/prey_marker_screen.md#all-host-depleted-libraries)
use [a separate manifest](../config/prey_marker_host_depleted.tsv) with all 151
Physalia and ten Nanomia septata libraries: 416,741,747 retained pairs. There is
no additional subsampling. The other 44 libraries are not part of these jobs.

| Job | State when submitted | Scope |
| --- | --- | --- |
| 13501446 | Running | COI reference-origin audit and candidate/control verification |
| 13501764 | Completed, exit 0 | Full fish cohort workflow on the nine existing synthetic pairs |
| 13501881 | Pending, dependency on 13501446 | 151 new COI screens; reuses ten Nanomia pilot searches; verifies all 161 |
| 13501886 | Pending, dependency on 13501446 | 156 new fish screens; reuses five Physalia pilot searches |
| 13501976 | Pending, dependency on 13501446 | Regenerate/check the COI cohort graph using the completed external references; no analyses |

Both expansion launchers inspect the scientific verification result and stop
before new read searches if it fails. Fish expansion also requires its completed
cohort control check. A successful scheduler exit for the audit is not sufficient
by itself. Scheduler state and the execution records, not this submission table,
establish completion.

The submitted source is frozen at
`data/results/prey_marker_submissions/2026-10-18-verified-expansion/`.
`submission.json` records job IDs and the SHA-256 of `submitted_source.tar.gz`.
Both jobs execute from that source directory. Its `data` and `logs` links point
to the existing project directories; later edits in the working tree do not
change these queued computations. Each job also preserves its executed source,
commands, versions, input identities, output hashes and exit status.

Verification before submission:

- All 16 focused COI tests passed; Ruff lint/format checks and launcher syntax
  checks passed.
- Fish cohort integration recovered the five expected candidate pairs from nine
  controls, exported both mates (ten reads), and completed competitive BLAST.
  Results and the passed check are under
  `data/results/fish_12s/cohort_validation/2026-10-18/`.
- All 130 saved fish pilot output hashes and the original COI pilot output hashes
  verified before reuse.
- Fish dry run: 156 new seed screens, 156 aligner comparisons, 156 paired unions
  and 156 competitive alignments. The five completed pilots remain inputs.
- COI scheduling dry run: 151 nucleotide screens, 302 translated mate screens,
  151 paired unions and nucleotide verifications, then 322 protein comparisons
  and 161 reports. The ten Nanomia libraries are not rescanned. This structural
  dry run used external-file placeholders under `/tmp` because the real curated
  panel and scientific gate were still being generated; no biological files or
  gate results were fabricated in the project.
- Verification and fish cohort graphs were generated and checked against their
  entry points. The COI cohort graph was initially generated with the same
  structural placeholders; job 13501976 regenerates/checks it against the real
  reference paths after the audit. These graph commands execute no analyses.

Neither workflow includes the original host-depletion, assembly, catalog,
mapping or manuscript rules. Fish filter membership and competitive alignments
will still require interpretation; a candidate is not automatically a fish
identification or evidence of predation. COI supported/ambiguous/unresolved
categories use the provisional criteria described in its guide.

One fish submission attempt also requested a dependency on completed control job
13501764; Slurm rejected that aged-out dependency. The successful submission uses
13501446 and checks the preserved passed-control result in its launcher.
