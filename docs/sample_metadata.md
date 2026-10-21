# Sample identifiers and collection metadata

The [publication sample table](../figures/cohort/sample_metadata.tsv) has one row
per analysis library. It can be updated without changing sequence-analysis
inputs. The frozen [manifest](../manifest.csv) and
[workflow sample table](../config/samples.tsv) retain the identifiers, read paths,
and settings used for the completed analyses.

## Identifier meanings

| Column | Meaning |
| --- | --- |
| `library_id` | Stable, unique library key, such as `Ahuja2024:CWD16`. Use this to join metadata additions. |
| `analysis_sample_id` | Existing workflow handle, such as `Ahuja2024__CWD16`. This is called `sample_id` in `config/samples.tsv` and detection tables. |
| `specimen_id` | Unchanged legacy identifier from the analysis manifest. Some values are museum vouchers; others are collection labels. |
| `voucher_id` | Museum catalog identifier including its institution and collection prefix, where known. YPM IZ records use `YPM:IZ:number`; other museums retain their source prefix. A field label alone is not a museum voucher. |
| `biosample` | NCBI BioSample accession. This identifies the submitted biological sample, separately from its SRA sequencing run. |
| `collection_id` | Field collecting-event identifier, such as an ROV dive and sampler label. Several tissue samples may share it. |
| `sample_id` | Source tissue/sample label, such as `CWD16`, `NA35`, or `IZT01143`. It is not the workflow handle and need not be globally unique. |
| `tissue_tube_id` | Additional dissection/tissue-tube label when a study records it separately. |
| `specimen_voucher` | Compatibility alias of `voucher_id`. |

For Church et al., source `YCGA_ID` and `diss_tube_ID` labels are supplemented
and corrected by Samuel Church's October 19 author review. His explicitly
confirmed CWD and S238 tissue labels take precedence; museum-derived submitted
labels retain the institution prefix. These labels are distinguished from a
separate tube number in the notes. Eight specimens never had a separate sample
label; their `sample_id` stays blank and their review decision is `not_assigned`.
Sequencing-center tracking codes and unconfirmed general specimen labels are
not substituted. Analysis identifiers remain unchanged.

## Collection and tissue fields

Latitude and longitude are decimal degrees, with south and west negative.
`location` includes a named collection site when available; `collection_site`
retains that site separately and `locality` retains the prior broad location.
Collection-event strings retain the source spelling, except the explicitly
recorded `Guam_104465`/`Guam-104465` spelling reconciliation for one voucher.

Depths are in metres. An explicit point depth takes precedence. When it is
missing and both minimum and maximum are available, `depth_m` is their arithmetic
mean and `collection_depth_basis` is `interval_midpoint`. Bounds and original
depth text remain available. A midpoint is a derived value, not a measured point
depth. A single bound is insufficient to calculate a midpoint. The earlier
[Physalia surface convention](metadata_depth_conventions.md) remains in force.

When a collection event includes animals sampled at multiple recorded depths,
an explicitly authorized representative mean is labeled `event_mean`.
Its minimum and maximum retain the observed event range. WS6 uses this basis;
its six source depths are retained in a
[dated table](../data/metadata/manual-review/2026-10-21_ws6_collection_depths.tsv).
This arithmetic mean is distinct from the midpoint of an interval and from an
individual specimen's measured depth.

`tissue` contains the source description. `tissue_status=recorded` means a value
is present; `not_recorded` means it still needs information. A study's general
extraction method is not substituted for a missing specimen-level tissue record.

`notes` contains human-facing explanations of missing identifiers or tissues,
derived depths, and known source disagreements. For example, NA19's note records
the decision to use 296 m despite another source reporting 297 m.
`analysis_notes` retains the older computational annotations, such as
population-cluster assignments, under a separate name. These are not evidence
about missing specimen metadata.

## Review and update

### Manual review worksheet

The [October 19 Church responses](../data/metadata/manual-review/church2025-2026-10-19/2026-10-19_church2025_review_response.md)
have been applied. They fill 24 previously missing values (16 tissues and eight
sample labels), resolve eight blank labels as never assigned, correct three
coordinate pairs, confirm 11 collection dates, and correct 22 existing sample
labels. All 70 original Church worksheet entries now have decisions. The
[application and remaining-work summary](../data/metadata/manual-review/church2025-2026-10-19/README.md)
separates 18 still-unresolved Church entries from eight never-assigned labels
and three retained dates with limited support. The counts in the table below
describe the original questions, not the current open workload.

On October 21, Casey Dunn resolved NA19 (V3698-D2) to **296 m**, accepting the
Ahuja2026 Table S8 point depth and treating the 1 m difference from the 2024
supplement/BioSample as immaterial. The decision is recorded in the
[Ahuja2026 worksheet](../data/metadata/manual-review/2026-10-17_ahuja2026_review.tsv)
and applied through the publication overlay. Casey subsequently supplied WS5's
D325-SS12 depth as **246 m** and directed use of the mean of six D326-SS1 depths
for WS6: **328.589662 m**, with observed bounds **295.8355–350.06003 m**.
The source sum is 1971.53797 m across six records; the stored mean is rounded to
six decimal places. The supplied log's `Nanomia bijuga` label is preserved in
the source table; the current host identification remains `Nanomia septata`.
The previous 2854.16 and 2563.33 values are superseded and retained in the
review history. All three Ahuja2026 depth questions are resolved. Earlier dated
question drafts and audit snapshots retain the original three-depth checklist.

The October 21 review also resolves seven entries in the
[Ahuja2024 worksheet](../data/metadata/manual-review/2026-10-17_ahuja2024_review.tsv):
CWD11 retains coordinates 40.31, -68.13; CWD131 retains date 2003-03-15;
CWD144's date is corrected to 2003-01-13; CWD16 retains depth 10 m; and NA8
retains depth 1299 m. CWD131's location is Gulf of California, Mexico, with
country Mexico; CWD146 retains date 2003-03-13. Notes retain the source
disagreements and record Casey Dunn's decisions. CWD16's source intervals
remain as provenance, with the confirmed point depth taking precedence.
The requested tissue, depth, date and location fields are now populated with
their flagged source disagreements resolved across the two Ahuja studies.
Separate voucher and collection-label questions remain open.

Use the study-specific worksheets below to supply corrections and decisions.
Each has an accompanying plain-text question draft that can be forwarded with
the worksheet to the study author. The drafts name every affected specimen,
give the existing and conflicting values, and describe the evidence needed.

| Study | Missing values | Catalog comparisons | Questions to forward | Worksheet to edit |
| --- | ---: | ---: | --- | --- |
| Church2025 | 50 | 20 | [Church questions](../data/metadata/manual-review/2026-10-17_church2025_questions.txt) | [70 review entries](../data/metadata/manual-review/2026-10-17_church2025_review.tsv) |
| Ahuja2024 | 11 | 13 | [Ahuja 2024 questions](../data/metadata/manual-review/2026-10-17_ahuja2024_questions.txt) | [24 review entries](../data/metadata/manual-review/2026-10-17_ahuja2024_review.tsv) |
| Ahuja2026 | 3 | 0 | [Ahuja 2026 questions](../data/metadata/manual-review/2026-10-17_ahuja2026_questions.txt) | [3 review entries](../data/metadata/manual-review/2026-10-17_ahuja2026_review.tsv) |

These are exact partitions of the [original October 17 review worksheet](../data/metadata/manual-review/2026-10-17_metadata_review.tsv):
64 missing-field entries and 33 flagged museum comparisons across 81 libraries,
prepared from the metadata at commit `cdd047c`. Keep that combined worksheet as
the initial checklist; enter responses in the study-specific files to avoid
maintaining two copies of an answer. The worksheets are human review documents,
not workflow inputs, and the exporter does not read them automatically.

CWD16 and NA19 occur in both Ahuja studies. Both question drafts cross-reference
them, but each has one response location: CWD16 in the Ahuja2024 worksheet and
NA19 in the Ahuja2026 worksheet. The Ahuja2026 draft also lists the 20 known YPM
vouchers whose public availability remained unconfirmed in the October 17 audit.
That catalog-release follow-up is separate from the three depth entries; public
URLs or dated status updates can be returned in prose.

Open it as a tab-separated table, importing identifiers as text. Edit these four
columns and retain the other columns as the record of what was reviewed:

| Column | What to enter |
| --- | --- |
| `decision` | `set_value` for an addition/correction; `keep_current` when the existing value is supported; `not_assigned` when an identifier was never assigned; `unresolved` when the available evidence cannot settle the question. Blank means not yet reviewed. |
| `reviewed_value` | The supported value for `set_value`. Otherwise leave blank. For coordinates, give both latitude and longitude with labels. For depth, give units and distinguish an actual depth from a reported range. |
| `evidence` | A record URL, notebook page, spreadsheet row, or named person's confirmation with date. For a group confirmation, specify which specimens it covers. |
| `review_notes` | A short explanation, including why a value remains unknown, why one source is preferred, or whether two event labels refer to the same event. |

Across the original three worksheets, 64 entries concern missing values.
After the reviews above, unresolved tissue and BioSample values remain for the
two OA specimens; no numeric depths are missing in the 205-library table.
The other gaps concern voucher, collection-event, and tissue/sample identifiers.
The other 33 entries compare existing study metadata with museum records.
A point outside a museum locality polygon is a discrepancy to review, not proof
that the study coordinates are wrong. An alternate event label
can be a valid synonym rather than an error.

Keep missing values blank rather than inserting guessed identifiers or the word
`unknown` as a biological value. A documented unresolved or never-assigned value
is a valid review outcome. Dates should be unambiguous (`YYYY-MM-DD` when a full
date is known), coordinates should use signed decimal degrees, and depth should
include units. If only a depth range is supported, report that range and let
curation calculate its midpoint. Do not average conflicting point depths.

After review, reconcile the decisions into the durable publication overlay,
preserve evidence for superseded values, update human-facing notes, and rebuild
the export as described below. In particular, the overlay's `original_value`
must refer to the frozen baseline rather than the worksheet's current value.
The generated publication table and the dated source-audit files are references,
not the place to enter review decisions.

### Durable updates and export

The [missing-field list](../figures/cohort/sample_metadata_missing.tsv) identifies
each unresolved value among voucher, BioSample, event, sample, coordinates,
depth, location, and tissue. Blank sample/event/voucher fields can represent
identifiers that were never assigned; do not invent them to make a row complete.
The [cell provenance](../figures/cohort/sample_metadata_provenance.tsv) gives
sources for those fields.

Durable publication edits belong in
[`data/metadata/publication_metadata_updates.tsv`](../data/metadata/publication_metadata_updates.tsv).
Each row records `library_id`, `field`, `original_value`, `value`, `source`,
`source_locator`, and `notes`. `original_value` is checked against the frozen
baseline (new fields start empty; `voucher_id` starts from `specimen_voucher`,
`location` from `locality`, and `analysis_notes` from the manifest's `notes`).
Human-facing `notes` starts empty; explanations of gaps and depth conventions are
then added automatically, alongside any sourced notes in the overlay.
Multiple agreeing sources are allowed;
conflicting values stop the export. When revising an already curated field,
reconcile its existing overlay rows rather than adding a conflicting value.
The older `sample_metadata_updates.tsv` is part of the frozen analysis baseline
and is not edited for this publication-only work.

Rebuild only the metadata with Python 3.10 or later; no additional packages,
sequence files, workflow invocation, or cluster jobs are required:

```bash
python3 scripts/export_sample_metadata.py
```

This writes the publication table, cell provenance, missing-field list, and
`sample_metadata_export.json`, which records input/output checksums. The cohort
figure builder uses the same exporter so a later figure rebuild preserves these
edits. Metadata export alone does not refresh existing plots, collection-context
summaries, or statistical results. Their provenance continues to identify the
inputs used when they were generated.

## October 2026 enrichment

The retained source snapshot is
[`data/sources/metadata_enrichment/2026-10-12-publication/`](../data/sources/metadata_enrichment/2026-10-12-publication/).
It includes the Church repository's specimen and sample tables at commit
`679ac669b4c44f64384373ca7adf15954f055231`, current public SRA experiment packages,
retrieval URLs/checksums, and the supplementary workbooks. Fresh downloads of the
Ahuja 2024 source workbook and [Ahuja 2026 Table S8](https://doi.org/10.1371/journal.pone.0351247.s008)
match the previously retained workbooks byte for byte. Their earlier reviewed
DNA-row matches supply the identifiers; RNA samples and out-of-cohort records
remain excluded. The additional Nanomia repository SRA workbook concerns genome
and other assays and supplied no new cohort values.

The initial enrichment's direct Peabody catalog requests returned an application
shell rather than specimen records. That enrichment used study sources. The
subsequent LUX catalog audit and metadata comparison are described below.

The update adds 18 BioSample/run accessions and 18 explicit tentacle records from
SRA. Coverage is 192/205 vouchers, 203/205 BioSamples, 194/205 collection IDs,
188/205 sample labels, 205/205 coordinate pairs and locations, 202/205 numeric
depths, and 187/205 tissue records. Fourteen depths are interval midpoints.
The remaining tissues are the 16 WAM specimens and two OA specimens; the two OA
specimens also lack BioSample accessions. NA19's conflicting depths and WS5/WS6's
unresolved depth meaning/units remain held, as documented in the
[previous source review](../data/metadata/enrichment/2026-10-04/README.md).

Reproduce the source proposal separately from the curated overlay:

```bash
python3 scripts/prepare_publication_metadata.py --output-dir /tmp/metadata-proposal
```

This uses retained source tables and XML with Python's standard library. The
Church S7 TSV was extracted from its retained workbook using openpyxl 3.1.5.
The accepted proposal, source-join audit, prior publication table, and checks
are retained in
[`data/metadata/enrichment/2026-10-12-publication/`](../data/metadata/enrichment/2026-10-12-publication/).

## October 17 LUX comparison

A fresh public-record audit again verified 148 of the 168 YPM vouchers in LUX.
The same 20 vouchers remained unconfirmed, including all the previously missing
NA28–NA38 and WS2–WS10 records. GBIF returned no exact matches. These are dated
search findings, not conclusions about internal museum access restrictions.
The [audit](../data/metadata/enrichment/2026-10-17-public-vouchers/summary.json)
retains requests, responses, and positive controls.

The [museum comparison](../data/metadata/enrichment/2026-10-17-lux/README.md)
reviews those 148 specimens and their 106 linked locality records. It adds 49
collection-site descriptions and corresponding location refinements, 26 explicit
country/territory entries, 22 collection methods/platforms, and depth bounds for
19 specimens. Existing point depths are retained; these additions create no new
midpoints. Museum voucher and matching sample/event labels gain corroborating
provenance. Optional fields, depth bounds, and human-facing notes are now included
in the exported cell-provenance table as well as the durable overlay.

Fourteen date disagreements, ten coordinate discrepancies, two depth
disagreements, six alternate field labels, and one inconsistent country/locality
combination are retained in a [review table](../data/metadata/enrichment/2026-10-17-lux/review_needed.tsv).
The existing study measurements and identifiers remain unchanged. Museum field
numbers and personal specimen numbers are not automatically reinterpreted as
missing event or tissue labels; museum locality numbers such as `IZS.34871` are
also kept out of `collection_id`. Locality polygons are checked for compatibility
with study coordinates, never converted to substitute point coordinates.

The original missing-field counts remain unchanged. The museum records do not
resolve NA19's conflicting depth or provide the still-missing tissue/BioSample
information. All frozen analysis inputs and their modification times are
preserved. Only the metadata export is rebuilt in this update.
