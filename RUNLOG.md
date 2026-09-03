# Run log

## 2026-09-03 — Phase 3 benchmark accepted; cohort-input gate ready

- Final benchmark recovery controller `10151500` completed all five remaining
  jobs in 3m02s. The 12-row decision table and independent benchmark validator
  report PASS for six libraries, three strategies, and all three host routes.
- CheckM2 assessed 39 bins. One strategy produced no bins, and six bins across
  three strategies were explicitly recorded as unassessable because they had
  no DIAMOND annotations; this state is distinct from a low-quality call.
- Locked cohort strategy: whole-genome reference depletion for `P_physalis`
  and `N_septata`; deterministic 25 M-pair trimmed-read input (seed 20260901)
  for reference-free libraries. Kraken-nominated assemblies had 88--100% host
  carryover in reference-bearing cases, while fixed effort recovered much more
  viral sequence in both reference-free comparisons.
- Locked objective assembly floor: at least 175,000 prepared read pairs. The
  lowest benchmarked input showing multi-domain recovery contained 176,589
  pairs. Input preparation is separated from assembly so the resulting
  eligibility table can be inspected and frozen before cohort MEGAHIT runs.
- The `phase3_inputs` stage implements the locked routing for all 205 libraries,
  aggregates retained-pair metrics, and validates exact cohort membership and
  eligibility. A three-study/three-route real-read fixture passed end to end
  and its clean dry run reports nothing to do.
- The production dry run contains exactly 412 jobs: 205 regeneration-only trim
  jobs, 205 cohort-input jobs, aggregation, and validation. Accepted Phase-2
  screens and the Phase-3 benchmark are not rerun. Next command after the
  implementation checkpoint is
  `bash scripts/submit_workflow.sh phase3_inputs config/config.yaml`.
- Controller `10151716` was submitted from implementation commit `946d922`.
  It is running on `week`; the first wave of 31 trim workers was accepted and
  started without an immediate workflow or scheduler failure.

## 2026-09-01 — Phase 3 assembly-strategy benchmark ready

- Phase 2 remains accepted: 205/205 screen validations, cohort sentinel PASS,
  independent validator PASS, and no active workflow jobs.
- The Phase-3 gate is frozen at six libraries and 12 paired strategy rows in
  `config/phase3_benchmark.tsv`. It spans high/low screen signal, all three
  studies, and all host routes. Reference-bearing libraries compare audited
  whole-genome depletion with Kraken-nominated reads; reference-free libraries
  compare Kraken-nominated reads with a deterministic 25 M-pair whole-read
  subsample (seed 20260901).
- Each strategy is evaluated by MEGAHIT assembly yield, barrnap rRNA markers,
  geNomad viral recovery, MetaBAT2 plus CheckM2 MAG recovery, host carryover
  where a conspecific reference exists, and measured runtime/memory.
- Static gate: manifest freeze PASS; independent Phase-2 cohort validation
  PASS; Phase-3 design/tool/database validation PASS; lint PASS; production
  dry run resolves 88 jobs. Runtime fixture tests pass for MEGAHIT, reference
  depletion, barrnap, geNomad, BWA back-mapping, MetaBAT2, and the no-bin
  CheckM2 path.
- Planned target: `phase3_benchmark` from the implementation checkpoint using
  `bash scripts/submit_workflow.sh phase3_benchmark config/config.yaml`.
- Controller `10089149` was submitted from checkpoint commit `818e747`; the
  persistent controller is running on `week` and began by materializing the
  pinned host-mapping environment needed for the two reference indexes.
- Live submission exposed a scheduler-policy boundary absent from dry-run:
  Yale rejects `week` workers requesting exactly 24 hours as too short for the
  partition. The immediately available MEGAHIT submissions were rejected;
  independent `day` workers were accepted. Correct the assembly and CheckM2
  week requests to 48 hours, let controller `10089149` and every accepted
  worker drain, then restart the identical target without overlapping jobs.
- Controller `10089149` exited 1 after 10h50m with 26/88 steps accepted; all
  workers drained. The corrected recovery dry run contained 62 jobs.
- Recovery controller `10105627` was submitted from correction commit
  `5739877`. All 12 MEGAHIT workers were accepted with 48-hour limits: 11
  started immediately and one initially waited only on `MaxCpuPerAccount`.
- This run produces a decision table only. Cohort assembly remains blocked
  until the benchmark result and objective inclusion threshold are accepted.

## 2026-09-01 — Cohort aggregation recovery checkpoint

- Recovery controller `10081702` completed all seven remaining library
  recoveries. All 205/205 libraries now have accepted screen-validation JSONs.
- Cohort aggregation worker `10085172` then reached its 2 GB memory allocation
  after 11m37s and was OOM-killed, so the controller exited 1 and the final
  `screen_cohort` sentinel was not created.
- Recovery resource decision: raise only `aggregate_phase2` to 8 GB through a
  controller-side scheduler override. Analytical parameters and accepted
  sample outputs are unchanged; `config/config.yaml` remains untouched.
- The recovery dry-run contains exactly two jobs: cohort aggregation and the
  final cohort sentinel.
- Controller `10087122` successfully published the 205-library aggregation
  with 2,282,231 nomination rows. Final validator `10087328` then exhausted
  its original 1 GB allocation while loading the 235 MB nomination table.
  Add an 8 GB scheduler-only override for `screen_cohort`; the next recovery
  dry-run should contain only that sentinel job.
- Final recovery controller `10087353` completed with exit 0. Sentinel worker
  `10087355` completed in 17 seconds with 2.44 GB peak RSS, and
  `data/results/stages/screen_cohort.done` reports PASS for 205 samples, three
  studies, three host routes, and 2,282,231 nomination rows.
- Independent `scripts/validate_phase2.py --scope cohort` reports PASS with
  zero errors. A final `screen_cohort` dry-run reports nothing to do.

## 2026-09-01 — Cohort resource-recovery checkpoint

- Replacement controller `10075623` resumed the 75-job remainder and completed
  56 jobs before exiting under `--keep-going`; 198/205 libraries now have
  accepted screen-validation JSONs and no controller or worker remains active.
- Full-cohort evidence exceeded the pilot envelope. PhyloFlash jobs for
  Ahuja 2026 NA35, NA38, and Ahuja 2024 NA10 reached the three-hour limit;
  Ahuja 2026 NA37, WS3, and NA31 reached the 24 GB memory limit. The uncapped
  Ahuja 2024 CWD11 trim also reached the 24 GB memory limit.
- Recovery resource decision: controller submissions override `trim_reads` to
  48 GB and `phyloflash_screen` to 48 GB / 6 h. These are scheduler-only
  corrections: analytical parameters and accepted outputs are unchanged, and
  `config/config.yaml` is deliberately not refreshed because doing so would
  invalidate the run snapshot for completed libraries.
- A clean recovery dry-run should contain 19 jobs: one trim, one Kraken/Bracken,
  one sylph, seven phyloFlash, seven sample validators, aggregation, and the
  final cohort sentinel.

## 2026-08-30 — Active cohort controller handoff checkpoint

- Active full-cohort controller: `9903780`, target `screen_cohort`, submitted
  from accepted gate commit `2f77618`; it is running on `day` with a hard
  one-day QOS limit and is scheduled to end at **2026-08-31 20:38:46 EDT**.
- Controller script commit `3f3b1f6` moves future controllers to `week` with a
  seven-day limit. This does not alter controller `9903780` or its workers.
- Do not cancel and overlap the current controller with a replacement. If it
  times out before `data/results/stages/screen_cohort.done` exists, let all
  submitted `sa.*` workers drain, then follow the controller handoff procedure
  in `WORKFLOW.md`. Only unlock after confirming that both the old controller
  and all of its workers are absent from `squeue`.

## 2026-08-30 — Phase 2 resource gate passed; cohort released

- Resource-gate controller `9744872`, commit `a71225d`, completed the refreshed
  43-job `screen_pilot` in 7h36m with exit 0. All eight trim jobs, 24 screens,
  eight content validators, aggregation, and the final sentinel completed; no
  worker failed, timed out, or exceeded memory.
- The refreshed sentinel reports `status PASS`, eight samples, three studies,
  three host routes, and 91,622 nomination rows. Independent
  `scripts/validate_phase2.py --scope pilot` reports eight validated samples
  and zero errors.
- A clean `screen_cohort` dry-run retains the eight completed pilot libraries
  and resolves 987 remaining jobs: 197 each of trim, Kraken2/Bracken, sylph,
  phyloFlash, and sample validation, followed by cohort aggregation and the
  final sentinel.
- Gate decision: the pilot-selected allocations and 200 M-pair cap are accepted
  for cohort screening. The full cohort may be submitted from this checkpoint.

## 2026-08-30 — Phase 2 pilot accepted; cohort resource gate

- Controller `9735088`, commit `29b7660`, completed the revised 200 M-pair
  `screen_pilot` target in 3h22m with exit 0. All eight libraries passed trim,
  Kraken2/Bracken, sylph, phyloFlash, content validation, and aggregation; the
  stage sentinel reports `status PASS`, three studies, three host routes, and
  91,623 nomination rows.
- The pinned Conda correction was effective: all eight Kraken2/Bracken workers
  completed without the former missing-Python failure. No pilot worker failed
  or exceeded its allocation.
- Observed resource maxima and selected cohort requests are recorded in
  `data/metadata/phase2_pilot_resources.tsv`. Kraken memory increases from 96
  to 104 GB (23% over the 84.2 GB peak). Other requests are reduced with
  margin because their pilot maxima were 15.8 GB RSS / 75 min for trim, 16.2
  GB / 20 min for sylph, and 16.6 GB / 81 min for phyloFlash.
- Gate decision: accept the 200 M-pair cap and proceed to the 205-library
  `screen_cohort` only after lint, frozen-manifest validation, pilot validation,
  capped-read fixture, and cohort dry-run pass from the resource-setting
  checkpoint.
- The resource edit refreshes the run snapshot by design, so the identical
  eight-library pilot will be rerun under the selected allocations before the
  cohort is released. This is the final resource gate, not a change to analysis
  parameters or the pilot sample set.

## 2026-08-29 — Phase 2 pilot revises compute ceiling to 200 M pairs

- Controller `9721858` completed all possible work in 7h59m: seven of eight
  pilot libraries passed trim, all three screens, and content validation.
  `Ahuja2024__CWD19` was OOM-killed during uncapped fastp at 384,583,011 raw
  pairs with 32 GB; the controller consequently exited 1 and did not aggregate.
- The 400 M ceiling reduced the planned cohort from 40.13 B to only 37.15 B
  pairs (7.4%) and left individual screen jobs processing approximately 380 M
  post-fastp pairs. It did not adequately control production cost.
- Decision: use a provisional **200 M-pair** production ceiling. This caps
  62/205 libraries, retains full depth for 143, reduces planned input to
  29.86 B pairs, and halves worst-case downstream screening effort. The
  100 M alternative would cap 163/205 libraries and is reserved for a
  pre-cohort sensitivity comparison because it may lose rare associate signal.
- Changing the locked config refreshes the run snapshot and invalidates prior
  cap-dependent trim/screen validation under the mtime-only restart policy.
  The revised eight-library pilot must pass before cohort submission.

## 2026-08-29 — Phase 2 full-library pilot ready

- Target: `screen_pilot` using `config/config.yaml`; eight libraries are frozen
  in `config/phase2_pilot.tsv` and span all three studies, all three host routes,
  capped/uncapped data, 31.5 M--1.184 B raw pairs, and one to five input pairs.
- Phase-2 screening is now independent of Phase-3 host handling. Per-library
  validators check synchronized trimmed counts, cap/provenance reconciliation,
  expected Kraken2/Bracken/sylph columns, phyloFlash archive contents, and
  screen provenance. Aggregation emits library QC and nomination-only evidence
  from Bracken, sylph, and phyloFlash.
- Scratch policy: trimmed production reads are temporary only after all three
  screens and validation succeed. Downstream priorities drain completed samples;
  exact capped sampling is restricted to one concurrent trim to bound temporary
  raw staging against 5.8 TB available scratch.
- Static acceptance: manifest freeze PASS; Snakemake lint PASS; pilot dry run 43
  jobs; cohort dry run 1,028 jobs; exact 1,000-pair cap fixture PASS; existing
  three-study fixture screen validation and 84-row nomination aggregation PASS;
  `scripts/validate_phase2.py --scope static` reports zero errors.
- Pilot controller `9721732` was submitted from commit `eb8bd30`. The first
  completed fastp run exposed an atomic-publication bug: reports created on
  `/vast` could not be renamed directly onto `/gpfs` (`EXDEV`). The data and
  fastp results themselves were valid. The doomed controller and remaining
  workers were canceled, and five enumerated incomplete scratch directories
  (409 GB) were removed; these contained no accepted output and are fully
  regenerable.
- Fix: atomic publication now falls back to a destination-filesystem temporary
  copy followed by rename. A direct `/tmp`-to-project cross-filesystem test
  passes. Next command after the fix checkpoint: restart the identical
  `screen_pilot` target; Snakemake will retain any accepted persistent output.

## 2026-08-29 — Phase 1 complete

- Workflow commit: `c702150` (`Build Phase 1 Snakemake workflow skeleton`).
- Target: `phase1_smoke` with `tests/config.fixture.yaml`; this is a 50,000-pair
  real-read fixture, not a cohort run.
- Coverage: Church 2025 / `P_physalis`, Ahuja 2024 / reference-free, and Ahuja
  2026 / `N_septata`.
- Static gates: Phase-0 frozen validation PASS; Snakemake lint PASS; production
  `screen_cohort` dry run resolves 1,234 jobs; SLURM-profile fixture dry run
  PASS; local fixture creation, fastp, and fixture-reference indexing PASS.
- Worker ceilings: Kraken2/Bracken 4 CPU / 96 GB / 2 h; sylph 4 CPU / 48 GB /
  2 h; phyloFlash 4 CPU / 16 GB / 2 h; fixture trim and host handling 4 CPU /
  8 GB / 1 h. Controller: 1 CPU / 4 GB / 24 h on `day`.
- Planned command:
  `bash scripts/submit_workflow.sh phase1_smoke tests/config.fixture.yaml`.
- Controller job `9721052`: FAILED before worker submission after 42 seconds.
  Cause: the generic-cluster helper ran under system Python and could not import
  `snakemake`; the one-CPU controller environment also exposed undesired worker
  thread scaling. Fix: run profile helpers with the workflow environment Python
  and give Snakemake an explicit scheduler-only core budget. Replacement
  controller: `9721122`, submitted from fix commit `285baf7`.
- Controller `9721122` worker finding: all three Kraken2 classifications ran,
  but jobs `9721143`, `9721146`, and `9721149` failed during atomic output
  handling because Kraken replaces `#` with `_1`/`_2`, producing a doubled
  underscore from the wrapper's original template. Corrected the template;
  locked Bracken to the database's available 150-mer distribution rather than
  an unsupported arbitrary post-trim mean length;
  completed upstream and independent screen outputs remain valid for restart.
- Restart controller `9721165`, submitted from wrapper-fix commit `797dad6`,
  completed in 1m59s. Worker jobs `9721167`--`9721169` (Kraken2/Bracken),
  `9721173`--`9721175` (host handling), `9721176`--`9721178` (sample
  validation), and `9721180` (stage validation) all completed with exit 0.
- The three validation JSONs report `PASS` and cover Church 2025 /
  `P_physalis`, Ahuja 2024 / reference-free Kraken-nominated reads, and Ahuja
  2026 / `N_septata`. The final sentinel reports three samples, three studies,
  and three host routes.
- Observed peak RSS informs the Phase-2 pilot: Kraken2/Bracken 84.3 GB, sylph
  14.8 GB, phyloFlash 4.7 GB, fastp <0.4 GB, and fixture host handling <0.24
  GB. These fixture values do not replace the required 6--10-library resource
  pilot on full libraries.
- Acceptance: `bash scripts/check_phase1.sh --require-smoke` PASS; Phase-0
  freeze PASS; Snakemake lint PASS; production dry run 1,234 jobs; completed
  fixture dry run reports nothing to do; `data/metadata/phase1_validation.txt`
  reports zero errors. Artifact: `tests/work/stages/phase1_smoke.done`.
- No cohort target was submitted. Next: begin Phase 2 by implementing and
  dry-running a predeclared 6--10-library `screen_pilot` target; do not submit
  `screen_cohort` until that pilot's resource gate is accepted.

## 2026-08-29 — Phase 0 complete

- Repository state: Phase-0 metadata built; no cohort analysis submitted or run.
- Analytical manifest: 205 unique libraries / 205 specimens.
- Primary-study rows: Church 2025 = 151, Ahuja 2024 = 33, Ahuja 2026 = 21.
- Raw inputs: 313 selected R1/R2 file pairs; sequencing batch parsed from each
  Illumina header as instrument/run/flowcell/lane.
- Provenance: 208 records. Church YPM-IZ-104465/Ahuja NA22 is one analytical
  library; CWD16 and NA19 are distinct Ahuja 2024 libraries each reused in the
  Ahuja 2026 study.
- Corrections: fixed two stale Church paths; excluded derived `combined_R1/R2`
  files; corrected five exact 2x read-count errors caused by counting combined
  files with their component lanes.
- Resources: 10 host/screen/taxonomy/QC database resources inventoried with
  versioned checksum anchors; both host-reference scaffold audits checked.
- Validation: `status PASS`, checksum freeze checked, zero errors. Report:
  `data/metadata/manifest_validation.txt`.
- Freeze: `data/metadata/manifest.freeze.sha256`.
- Scale estimate: 40.13 billion raw pairs in 6.05 TB compressed; 12 libraries
  exceed the 400 M-pair ceiling. Estimated capped trimmed-read reservation with
  20% margin is 6.68 TB, larger than the 5.8 TB scratch space observed during
  planning. Phase 1 must stage libraries or reduce simultaneous scratch
  residency; this does not change the locked 400 M-pair ceiling.
- SLURM jobs: none.
- Next: Phase 1 — build and smoke-test the Snakemake workflow and SLURM profile.

Reproduce and verify Phase 0:

```bash
python3 scripts/build_manifest.py
python3 scripts/inventory_phase0_resources.py
python3 scripts/summarize_phase0.py
python3 scripts/validate_manifest.py
python3 scripts/freeze_manifest.py
python3 scripts/validate_manifest.py --freeze data/metadata/manifest.freeze.sha256
```
