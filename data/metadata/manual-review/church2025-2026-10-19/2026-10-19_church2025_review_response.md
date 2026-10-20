# Church2025 specimen metadata: review response

Responding to: 2026-10-17_church2025_questions.txt
Response date: 2026-10-19
Answers by: Samuel Church
Detailed answers: `2026-10-19_church2025_review_responses.tsv` (70 review entries plus 22 additional label corrections)
Main evidence source: `Physalia_2022_2023_Collecting_POPGEN.xlsx` (project collecting, dissecting, sequencing and submission records). Row numbers in the TSV refer to this workbook.

## Identifier policy

For every Physalia specimen, we use the museum voucher code, including the institution prefix (YPM-IZ-, WAM-, TMAG-, FM-, NIWA-), as the primary identifier. We deliberately avoid relying on tube or sample codes where possible. Many specimens never had a separate sample number. Where no separate tube label exists, the tissue/sample label is recorded as `not_assigned`, not filled in with the museum number or a sequencing-center tracking ID.

## Source priority

For dates and coordinates, we prioritized the original collecting record in the project spreadsheet over the museum catalog (LUX).

## Summary by question

1. **Tissue used.** All 16 WAM specimens: tentacle (subsamples provided alongside the whole specimens). OA-190531-1 and OA-Bluebottle: unknown.
2. **OA specimens.** Voucher, BioSample and tissue are unknown for both. The specimens were provided by Ozren and Alex. GON10891 is kept as OA-Bluebottle's field event; whether it is an event or a tissue identifier is unconfirmed.
3. **Field event IDs.** All 11 are unknown. FM, TMAG and WAM specimens were collected by others; event details may be available from those institutions' catalogs. The four YPM specimens have no event ID in any project record (collectors are named in the TSV).
4. **Tissue/sample labels.**
   - CWD numbers: YPM-IZ-111207 to 111210 (CWD_NZ31–34) and NIWA-173304/173305 (CWD_NZ35/36; formerly YPM-IZ-111211/111212).
   - Event names: YPM-IZ-111839/111840 (S238-015/S238-013).
   - None assigned: YPM-IZ-111156 to 111160 and 111163 (A. Schaeffer collections) and YPM-IZ-110434/110436 (L.M. Rojas collections, sent straight to the museum).
   - OA-Bluebottle: unknown.
5. **Collection dates.** All 11 current dates are kept, as they match the collecting record. For YPM-IZ-111006 to 111010, the catalog date (2023-11-08) is the DNA extraction date. Support is weaker for YPM-IZ-110460 (collector not confirmed) and YPM-IZ-111241/111243 (date appears only in the submission sheet).
6. **Coordinates.**
   - Six current values are kept.
   - Three are corrected to the raw collecting record: YPM-IZ-110613, 110906 and 110907. The current values do not match the raw record; for 110906 and 110907 this is a degrees–decimal-minutes conversion error. All three corrected values fall within the museum's mapped extent.
   - YPM-IZ-110631 and 110667: the original records confirm west (negative) longitude, so the museum geometry has a sign error.
   - YPM-IZ-110906 and 110907 are separate collections on different days, not positions along one trawl.

## Additional corrections (outside the 70 entries)

- TMAG tissue labels KM5632–KM5636 are incorrect and should not be reproduced. Use TMAG-K5632 to TMAG-K5636, which match the SRA/BioSample names.
- WAM tissue labels should always include the institution code (for example, WAM-Z88480).
- FM-16644: the museum number is the only identifier.

## Follow-up outside this review

- Update the NCBI BioSample coordinates for YPM-IZ-110613, 110906 (SAMN40621373) and 110907 (SAMN40621374).
- Add the tissue (tentacle) to the WAM BioSamples, where it is currently "missing".
- Report the longitude-sign error (110631, 110667) and the extraction-date-as-collection-date issue (111006 to 111010) to the YPM collection manager.
