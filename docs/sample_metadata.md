# Sample identifiers and collection metadata

The [publication sample table](../figures/cohort/sample_metadata.tsv) has one row
per analysis library (213 across four studies). It can be updated without changing sequence-analysis
inputs. The frozen [manifest](../manifest.csv) and
[workflow sample table](../config/samples.tsv) retain the identifiers, read paths,
and settings used for the completed analyses.

The expanded selection is [post_expansion_samples.tsv](../config/post_expansion_samples.tsv).
The exporter combines the original 205-library baseline and its reviewed overlay
with [Manko2026 specimen metadata](../data/metadata/manko2026/sample_metadata.tsv)
and [study manifest](../manko2026_manifest.csv). It validates exact expanded
membership; adding metadata does not modify either study's frozen sequence inputs.
For Manko2026, publication-overlay original values are the sourced specimen values;
original study identifiers remain protected.

## Nanomia identifications

The publication table distinguishes *Nanomia bijuga*, *N. cara*, *N. septata*,
and *Nanomia* sp. 1. Voucher-level identifications in the mitochondrial trees
of Ahuja et al. (2026; Fig 4 and S8 Fig), joined through their S8 Table, resolve
six specimens previously pooled as `Nanomia sp.`. The
[reviewed identifications and evidence](../data/metadata/manual-review/2026-10-23_nanomia/README.md)
record each assignment. Current names are applied through
[the publication host-name table](../config/publication_host_names.tsv);
original sequence-analysis manifests retain their deposited labels.

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

## Coverage and unresolved records

The [export record](../figures/cohort/sample_metadata_export.json) reports 213
libraries with 197 museum vouchers, 211 BioSamples, 202 collection identifiers,
196 tissue-sample labels and 211 tissue descriptions. Coordinates, numeric depth,
location and collection date are populated for all 213; populated values can
still have documented source conflicts. The
[missing-field table](../figures/cohort/sample_metadata_missing.tsv) identifies
remaining gaps, and [cell provenance](../figures/cohort/sample_metadata_provenance.tsv)
records their sources.

OA-190531-1 and OA-Bluebottle still lack verified public sequencing accessions,
BioSamples and specimen-level tissues. No general extraction method is substituted
for a missing tissue. Some other blank voucher, collection-event and sample labels
were never assigned. Eight Church specimens explicitly never had separate sample
labels, so those fields remain blank. Unknown or never-assigned identifiers must
not be invented to complete a row.

Reviewed depth values include NA19 at 296 m, WS5 at 246 m and WS6 at the event
mean of 328.589662 m. All earlier numeric-depth gaps are resolved. The
[depth guide](metadata_depth_conventions.md) explains the 150 assigned surface
values, 40 source records, 14 interval midpoints, one event mean and eight reported
point depths. Five Mańko et al. NCBI source records disagree with the published
collecting-event metadata; the [study guide](bargmannia.md#collection-metadata-source-conflicts)
records the retained values and their provenance.

## Review evidence

The dated study worksheets and source snapshots preserve decisions and original
questions. Their original question counts are not current missing-value counts.
Use the generated missing-field table above for the current gaps.

- [Church author responses](../data/metadata/manual-review/church2025-2026-10-19/README.md)
  document tissue, sample-label, coordinate and date decisions.
- [Ahuja 2024 worksheet](../data/metadata/manual-review/2026-10-17_ahuja2024_review.tsv)
  records reviewed dates, locations, coordinates and depths alongside unresolved
  identifier questions.
- [Ahuja 2026 worksheet](../data/metadata/manual-review/2026-10-17_ahuja2026_review.tsv)
  and the [WS6 source depths](../data/metadata/manual-review/2026-10-21_ws6_collection_depths.tsv)
  record the resolved depth questions.
- [Nanomia identifications](../data/metadata/manual-review/2026-10-23_nanomia/README.md)
  connect current host names to published trees and vouchers.
- [LUX museum comparison](../data/metadata/enrichment/2026-10-17-lux/README.md)
  retains the catalog evidence and source disagreements. Its unresolved counts
  describe that dated audit, before the author reviews.

A point outside a museum locality polygon is a discrepancy to assess, not proof
that the study coordinates are wrong. Museum field numbers and personal specimen
numbers are not automatically substituted for collection-event or tissue labels.
Keep dates unambiguous, use signed decimal-degree coordinates, and retain source
units and depth intervals. Do not average conflicting point depths.

## Durable updates and export

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
