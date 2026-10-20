# Church author responses applied October 20, 2026

The two October 19 responses from Samuel Church are retained unchanged here.
All 70 requested decisions and 22 additional label corrections were matched
by library and field. The durable publication overlay and the original Church
review worksheet now incorporate those decisions. `applied_changes.tsv` records
the changes; `metadata_before.tsv` preserves the publication table at checkpoint
`9d3ed92`, and `superseded_publication_updates.tsv` preserves replaced evidence.

## Gaps closed and corrections

- Tissue recorded as tentacle for all 16 WAM specimens; only the two OA tissues
  remain unknown in Church2025.
- Eight sample labels supplied: CWD_NZ31–34 for YPM-IZ-111207–111210,
  CWD_NZ35/36 for NIWA-173304/173305, and S238-015/S238-013 for
  YPM-IZ-111839/111840.
- Eight labels confirmed never assigned: YPM-IZ-110434, 110436, 111156–111160,
  and 111163. Their values remain blank, with explanations and `not_assigned`
  decisions. They do not need another search for a missing tube label.
- All 11 date comparisons received keep-current decisions; for 111006–111010,
  the museum date was the DNA extraction date. Three dates still have weaker
  support, listed below.
- Six coordinate pairs confirmed; three corrected:

| Specimen | Latitude | Longitude |
| --- | ---: | ---: |
| YPM-IZ-110613 | -30.884074 | 153.066594 |
| YPM-IZ-110906 | 26.466383 | -164.766667 |
| YPM-IZ-110907 | 26.004800 | -176.263200 |

The two latter specimens were collected on different days (August 23 and 22,
2023), not along one trawl. The existing negative longitudes for 110631 and
110667 are confirmed; the museum geometries have sign errors.

- Twenty-two existing sample labels corrected to include the institution:
  five TMAG, 16 WAM and FM-16644. These are submitted labels based on museum
  identifiers, not newly inferred tissue tube identifiers.

## Remaining work

Church's 50 original missing-field entries are now 24 filled, eight confirmed
not assigned, and **18 unresolved**:

- OA-190531-1 and OA-Bluebottle: tissue, voucher and BioSample (six entries).
- OA-Bluebottle: sample label (one). GON10891 remains in `collection_id`, but
  whether it describes an event or a tissue sample remains unconfirmed.
- Collecting-event identifiers for FM-16644, TMAG-K5632–K5636, WAM-Z88480,
  YPM-IZ-110694, 110906, 110907 and 111163 (11). The last four have no event ID
  in the reviewed project records; external museums may have additional records
  for their specimens.

Three dates are retained with limited evidence: YPM-IZ-110460 (collector
confirmation pending), and YPM-IZ-111241/111243 (submission-sheet evidence only).
The 24 Ahuja2024 questions and three Ahuja2026 questions are unchanged. They
include 11 missing Ahuja2024 vouchers, 13 catalog comparisons, and the depths of
NA19, WS5 and WS6. See [remaining_review.tsv](remaining_review.tsv) for all
remaining entries, including the three qualified Church dates.

The generated missing-field table still lists blank biological values, including
the eight never-assigned labels. It therefore contains 26 Church blanks; this
is not a count of unresolved review questions.

External follow-up remains: correct three NCBI coordinate records, add WAM
BioSample tissue descriptions, and report the museum longitude/date errors.
No external records or messages were changed by this application.

## Reproduction and scope

From checkpoint `9d3ed92`, run `python3 scripts/apply_church_metadata_review.py`
then `python3 scripts/export_sample_metadata.py`. The application deliberately
fails if response baseline values differ; do not reapply it to already corrected
metadata. The exporter alone can be rerun safely.

All 205 library memberships, `analysis_sample_id`, `specimen_id`, read paths and
analysis manifests are unchanged. This metadata-only update does not rerun
sequence analyses. Existing geographic figures/statistics retain their previous
input provenance and should be refreshed when the expanded-cohort reporting is
rebuilt. The separate Bargmannia metadata are outside these Church responses.
