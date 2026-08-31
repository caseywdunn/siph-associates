# Run log

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
