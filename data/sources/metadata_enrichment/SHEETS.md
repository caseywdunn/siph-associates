# Primary-sheet metadata extraction

The user supplied two primary supplementary workbooks on 2026-10-04.
Original XLSX files are preserved in `data/sources/` under the existing binary
source ignore rule. `sheet_sources.json` records their checksums and column
definitions. `sheet_raw_rows.json` preserves source rows and cell addresses.

Run from the analysis repository, using Python with openpyxl 3.1.5:

```bash
python scripts/extract_sample_sheet_metadata.py \
  --input-dir ../manuscript_siph_associates/tmp --repo .
```

The extraction writes candidates only; it never changes `manifest.csv`.
`sheet_metadata.tsv` uses the shared long-form merge schema. Its provenance
includes source file, sheet/cell address, source sample and matching IDs.
`sheet_metadata_pending.tsv` holds WS5/WS6 numeric depths for external review;
do not apply it automatically. `sheet_row_matches.tsv` records exclusions.

There are 78 source rows: 57 matched DNA rows covering 55 distinct cohort
libraries, 20 excluded RNA rows, and one DNA library outside this cohort
(SRR24493087/YPM-IZ-111756). Study aliases preserve the NA22/Physalia,
NA19 and CWD16 reuse mappings without duplicating libraries. An ambiguity
between identifiers excludes a row rather than choosing one identifier.

Keep depth intervals as minimum/maximum bounds. CWD16's 10 m record is
compatible with its 0–20 m interval. NA19 differs between primary sources:
297 m versus 296 m; neither source is silently replaced. The 2026 column
has explicit metre units in interval cells, but not in numeric cells;
normalization records this interpretation in `notes`. Raw depths are retained.
No taxonomy, sequencing metadata or detection outputs are changed.

Validation covered complete cohort joins, source checksums, normalized
coordinate bounds, all interval endpoints, the NA22/CWD16/NA19 alias cases,
and exclusion of RNA/HiFi records. Ruff formatting and lint checks passed.
