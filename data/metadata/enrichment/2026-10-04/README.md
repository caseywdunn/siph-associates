# Collection metadata enrichment, 2026-10-04

This is the source-enrichment snapshot before the subsequent authorized
Physalia surface-depth convention. Its counts and input checksums describe
that stage. The current overlay additionally includes the 150 assigned zero
depths documented in `docs/metadata_depth_conventions.md`; the source-enrichment
overlay and manifest are preserved alongside this README.

The user supplied two specimen spreadsheets and authorized filling missing
collection metadata using them and public BioSample/YPM records. The reviewed
update fills 895 cells in 185 of the 205 libraries. Every pre-existing
nonmissing value, library identifier, host assignment, input-read path, read
count, inclusion decision and sequencing batch is preserved.

## Coverage

| Metadata | Before | After |
| --- | ---: | ---: |
| Coordinate pairs | 172 | 205 |
| Numeric collection depth or interval | 0 | 52 |
| BioSample accession | 0 | 185 |
| Tissue description | 0 | 169 |
| Country/territory | 0 | 126 |
| Life stage | 0 | 5 |
| Named study ocean region | 172 | 173 |

Depths comprise 38 points and 15 intervals, with CWD16 represented in both;
there are 52 distinct specimens, including 51 of 54 non-Physalia specimens.
Point depths span 0–1349 m. The sole Physalia depth is NA22/YPM-IZ-104465,
0 m. The supplied sheets and retrieved BioSamples do not identify which other
Physalia were beach-collected. The user's zero-depth convention for documented
beach collections was therefore not applied to all Physalia by inference.
`depth_source` remains sequencing-read-count provenance;
`collection_depth_source` records collection-depth evidence separately.

`coverage_summary.tsv` reports each field; `changed_cells.tsv` lists each
addition. `metadata_merge_audit.tsv` retains agreement, fills and conflicts,
including source files/URLs, cells/attributes and original values. The accepted
overlay is `../../sample_metadata_updates.tsv`; rebuilding the manifest applies
it after library deduplication. It reproduces the promoted manifest exactly.

## Evidence and review decisions

The source directory is `data/sources/metadata_enrichment/` relative to the
analysis repository. Its spreadsheet extraction covers 57 DNA rows representing
55 distinct libraries. Twenty RNA rows and one out-of-cohort DNA row were
excluded. Original workbooks, raw rows, cell locators, checksums and alias joins
are retained. BioSample requests retrieved all 185 mapped accessions; 20
libraries have no public run accession in the current manifest. XML responses,
request URLs, timestamps, original attributes and exact accession checks are
retained under `external/`.

The GBIF check searched the Yale Peabody Invertebrate Zoology dataset
`854e35e6-f762-11e1-a439-00145eb45e9a` for all 168 complete YPM:IZ vouchers,
using both padded and unpadded catalog numbers. It returned no exact cohort
matches; positive controls confirmed the dataset and repeated-parameter query
worked. A complete 80-record Physalia search also had no cohort voucher matches.
This bounded query supplies no additions and does not establish that museum
records are unavailable elsewhere. GBIF requests and coverage remain under
`external/gbif/`.

Twenty-four explicit decisions normalize 22 YPM voucher spellings and exclude
two lower-precision CWD16 longitude candidates. The retained 2026 longitude,
−74.016667, rounds to the older −74.02 and fills a previously blank field.
No coordinates were averaged. The original values remain in
`review_decisions.tsv` and `reviewed/review_decision_audit.tsv`.

Unresolved issues remain visible:

- NA19: 297 m in the 2024 sheet/BioSample versus 296 in the 2026 sheet.
  Numeric depth remains blank; these values were not converted into a spurious
  collection interval.
- WS5 and WS6: 2854.16 and 2563.33 in the sheet and BioSamples. Their meaning
  and units need confirmation; raw values are retained but numeric depth is
  withheld. The BioSample submissions repeat the sheet rather than providing
  independent measurements. WS7–WS10 numeric depths use the spreadsheet's
  metre-labelled column context, recorded in their source notes; BioSample
  numbers without units were not independently converted.
- NA34–NA36: BioSample collection dates are 2011-06-16, versus 2024-06-16 in
  the sheet and existing metadata. The existing 2024 dates are preserved.
- Guam NA22/YPM-IZ-104465: the sheet and BioSample give west longitude,
  inconsistent with Guam and the existing east coordinate. Existing coordinates
  are preserved; the source discrepancy remains in the audit.
- CWD16's 10 m point and 0–20 m range are compatible and both retained.
  Their different raw strings are preserved in the audit, rather than choosing
  one original string. Other audit differences chiefly concern precision,
  locality granularity or formatting; existing values remain unchanged.

## Reproduction and validation

The original manifest and freeze are preserved in
`data/metadata/history/2026-10-04_before_sample_enrichment/`. To reproduce the
review, run from the analysis root, with new output directories:

```bash
python scripts/extract_sample_sheet_metadata.py \
  --input-dir ../manuscript_siph_associates/tmp --repo .
python3 data/sources/metadata_enrichment/external/prepare_biosample_requests.py
# Fetch the recorded URLs; existing raw responses permit offline normalization.
python3 data/sources/metadata_enrichment/external/normalize_biosamples.py
python3 scripts/review_sample_metadata.py \
  --candidates data/sources/metadata_enrichment/sheet_metadata.tsv \
    data/sources/metadata_enrichment/external/biosample_candidates.tsv \
  --decisions data/sources/metadata_enrichment/review_decisions.tsv \
  --output-dir /tmp/new-siph-metadata-review
python3 scripts/merge_sample_metadata.py \
  --manifest data/metadata/history/2026-10-04_before_sample_enrichment/manifest.csv \
  --candidates /tmp/new-siph-metadata-review/reviewed_candidates.tsv \
  --updates /tmp/nonexistent-siph-metadata-overlay.tsv \
  --output-dir /tmp/new-siph-metadata-proposal
```

Extraction uses openpyxl 3.1.5. Other metadata scripts use the Python standard
library. The review/merge scripts write proposals only; promotion copied the
verified proposed manifest and overlay to their canonical paths. Twenty-two
targeted tests cover exact identifiers, preservation, conflicting depths,
coordinates, source provenance and replay. An independent review confirmed the
merge. The existing manifest validator passed with 205 specimens, 313 paired
FASTQ records, 208 provenance rows, ten resources, two host-reference audits
and zero errors. `promotion_validation.json` records hashes and invariants.

The sole study-region addition is CWD16: blank to NW Atlantic. This changes its
descriptive label but not the permutation groups: it remains a one-library,
one-flowcell Nanomia bijuga stratum. Permutable membership remains 11 strata /
115 libraries overall and two strata / six libraries for Nanomia. Physalia
regional-model inputs are unchanged. Historical inference tables retain the
metadata used for their execution; sequence analyses and statistics were not
rerun. All 205 libraries remain descriptive; 136 are eligible for flowcell
inference, including 122 Physalia libraries.

Cohort publication exports are refreshed with:

```bash
MPLCONFIGDIR=/tmp/siph-associates-matplotlib \
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python figures/build_cohort.py
```

The current annotation freeze includes the accepted overlay and enrichment
sources/scripts. The historical executed-analysis freezes are preserved.
Source base revision is `0a7c3e85a0e8f4856a16e5fc75e6a0ec73ba7b19`, plus
uncommitted working-tree changes; checksums identify the actual artifacts.
