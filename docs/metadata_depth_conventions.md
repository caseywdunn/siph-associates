# Collection-depth conventions

All 213 specimens have a numeric depth in the
[current publication metadata](../figures/cohort/sample_metadata.tsv). The
`collection_depth_basis` field distinguishes measured or reported values from
derived ones. Source text and available interval endpoints remain in the table.

| Basis | Specimens | Meaning |
| --- | ---: | --- |
| `curator_assigned_surface` | 150 | *Physalia* with no recorded depth assigned 0 m |
| `source_record` | 40 | Depth retained from source metadata, including one *Physalia* at 0 m |
| `interval_midpoint` | 14 | Mean of recorded minimum and maximum when no point depth exists |
| `event_mean` | 1 | Explicitly reviewed mean of six collection depths for the WS6 event |
| `reported_point` | 8 | Point depths reported for the Mańko et al. specimens |

## Assigned surface values

For *Physalia*, a value of 0 m is assigned only if no point, interval or original
depth record is present. Explicit records take precedence. The convention gives
all 151 *Physalia* a zero depth, but it does not establish that a specimen was
beach-collected or separately measured at the surface. Original depth text is
not invented for assigned records.

The original proposal and source audit are retained in
`data/sources/metadata_enrichment/physalia_depth_scoring.tsv` and
`data/metadata/enrichment/2026-10-04-physalia-surface/`.

## Intervals and event means

A recorded point takes precedence over an interval midpoint. If the point is
missing and both bounds exist, the exporter uses their arithmetic mean and
labels it `interval_midpoint`. A single bound is insufficient. All bounds remain
available, and the full interval determines the collection-context depth stratum.

NA19 uses the reviewed 296 m point and WS5 uses 246 m. WS6 instead uses
328.589662 m, the arithmetic mean of six depths for collecting event D326-SS1,
with observed bounds 295.8355–350.06003 m. The
[six source depths](../data/metadata/manual-review/2026-10-21_ws6_collection_depths.tsv)
sum to 1971.53797 m. This is an event mean, not an interval midpoint or a measured
depth for the individual sequenced specimen. Collection-context tables label it
`range_with_derived_mean` and retain the full range.

## Interpretation and updates

Depth summaries retain the source basis. Comparisons of surface *Physalia* with
other hosts combine habitat, host identity and sampling design. Within
*N. septata*, locality and depth overlap. These data do not isolate depth effects.

Use the [metadata guide](sample_metadata.md#durable-updates-and-export) to review
or update depths. A metadata export alone does not regenerate downstream
summaries or statistical results; the [reporting workflows](expanded_reporting.md)
consume the reviewed values. Earlier source snapshots retain their recorded
values and provenance.
