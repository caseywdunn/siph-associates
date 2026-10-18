# Crustacean COI pilot execution

Work is isolated on `fish-12s-screen`. The current biological scope, inputs and
expansion policy are in [the prey-marker plan](../../docs/prey_marker_screen.md);
the executable method is documented in [the COI guide](../../docs/crustacean_coi_pilot.md).

## Revised run

On 2026-10-17 Casey requested all retained host-depleted pairs for reference-bearing
hosts, with initial fish 12S screening in *Physalia* and COI in *N. septata*.
`config/crustacean_coi_pilot.tsv` now contains all ten *N. septata* libraries,
2,152,939 pairs in total. No new read filtering, host mapping or subsampling occurs.
Job **13499718** completed references, synthetic controls and these inputs; real reads
depended on the control gate passing. Execution evidence is retained in
`data/results/crustacean_coi/executions/13499718/` and per-library `summary.json`.

The reference build and all twelve end-to-end paired controls passed.
[Control results](control_results.json)
show that host, fish and insect COI can be recruited by both methods but are
correctly assigned as non-target during competitive alignment. The synonymous
copepod control is recovered by both methods; it does not demonstrate a
translated-only sensitivity gain. Random and poly-G controls are not recruited.

Code checks passed before submission: Ruff formatting/lint, nine focused pytest
checks of genetic-code translation, COI annotation and competitive assignment,
and a dry run of the isolated DAG. Existing analyses are not dependencies.

Reference extraction retained 17,807 annotated RefSeq COI records and 128,274
MetaZooGene barcode records before frame validation. It excluded 28 features
whose RefSeq records omit the actual sequence, 314 annotated features with
unsuitable length or ambiguity, and the per-group barcode exclusions recorded
in `panel/2026-10-17/extraction.json`. Final coverage follows frame validation.
The validated panel has 108,778 target records and 36,459 competitors;
[coverage](reference_coverage.tsv) includes 919 copepod species labels and 68
euphausiacean labels. These source labels are not an independently reconciled
species count. Frame validation excludes 768 unsupported barcode frames and
76 translations with stops or insufficient length. Compact reference summaries
and the original [download receipt](reference_downloads.json) are retained here.
Download identities and hashes are recorded in the reference receipt. The
RefSeq release listing dates the downloaded GBFF to 2026-09-03; the downloaded
MetaZooGene FASTA headers say `MZGdb_v2025-m10-01`.

## Superseded attempts

The [earlier manifest](superseded_raw_pilot.tsv) selected the first ten million
pairs from one raw lane each of NA4 (*Praya*) and NA37 (*Nanomia* sp.), based on
existing copepod SSU evidence. No real-read screening completed for that design.

- 13499572: extraction failed on a RefSeq feature with undefined nucleotide
  content; the parser now counts and skips such records.
- 13499579: extraction completed; the broad barcode-frame search was stopped
  after about five minutes to reduce redundant comparisons. Frame validation
  now uses unique annotated crustacean COI proteins; competitive identification
  retains the full non-target panel.
- 13499614: failed on the cancelled job's stale lock; scheduler completion and
  COI-only lock contents were verified before unlocking.
- 13499669: the old two-library scope was cancelled when Casey revised the plan.
  No input subset had been created. Source snapshots and logs remain in the
  corresponding execution directories.

The main workflow, specimen identifiers and accepted manuscript results were
not changed. Raw full-depth reference-free reads will be considered only after
the revised pilot works and resource use is assessed.

## Completed pilot and follow-up

Job 13499718 completed in 15:04 with exit 0, screening 2,152,939 pairs and
recruiting 92 candidate pairs. Its original 28 non-target, 62 unresolved, one
ambiguous and one provisional crustacean labels are preserved. Subsequent
[validation](validation.md) identified a ranking artifact and questionable
reference labels; those initial categories are not accepted prey detections.
