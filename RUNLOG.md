# Run log

## 2026-08-29 — Phase 1 smoke planned

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
  controller: pending validation and submission.
- Acceptance artifact: `tests/work/stages/phase1_smoke.done`, followed by
  `bash scripts/check_phase1.sh --require-smoke`.

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
