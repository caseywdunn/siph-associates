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

For Church et al., `sample_id` uses the source `YCGA_ID`, or `diss_tube_ID` when
the sequencing-center label is absent. Both identify the submitted tissue/sample;
the original column is retained in the provenance. General `specimen_ID` labels
are not silently treated as tissue identifiers. Repeated sequencing entries are
accepted only when their labels agree. The published NIWA and TMAG museum-ID
corrections are applied when matching source records.

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

`tissue` contains the source description. `tissue_status=recorded` means a value
is present; `not_recorded` means it still needs information. A study's general
extraction method is not substituted for a missing specimen-level tissue record.

`notes` contains human-facing explanations of missing identifiers or tissues,
derived depths, and known source disagreements. For example, NA19's depth is
blank because two sources report 296 and 297 m. `analysis_notes` retains the
older computational annotations, such as population-cluster assignments, under
a separate name. These are not evidence about missing specimen metadata.

## Review and update

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
