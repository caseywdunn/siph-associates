# Run log

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
