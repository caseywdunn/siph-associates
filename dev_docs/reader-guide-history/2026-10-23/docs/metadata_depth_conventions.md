# Collection-depth scoring convention

On 2026-10-04 the user instructed: "Can score the depth for all physalia
without explicit record as 0." This extends the preceding specimen-metadata
enrichment; it does not revise the source sheets or BioSample records.

Score `depth_m=0` for a Physalia library only when no point, interval or original
depth record is present. Preserve explicit depth records. The new
`collection_depth_basis` column distinguishes `curator_assigned_surface` from
`source_record`. A scored zero does not establish that a particular specimen
was beach-collected, and it is not presented as a separately measured depth.
Raw `depth_original` fields are not invented for scored records.

The convention adds 150 zero depths to the one explicitly recorded Physalia
zero. All 151 Physalia therefore have depth 0 m; 202 of 205 cohort libraries
have a point or interval. NA19, WS5 and WS6 remain unresolved. Point depths
and interval endpoints of all non-Physalia libraries are unchanged.

`scripts/score_physalia_depths.py` generates candidates from the preceding
manifest, and `scripts/merge_sample_metadata.py` applies them through the
existing reviewed overlay. Inputs and the proposal/audit are retained in
`data/sources/metadata_enrichment/physalia_depth_scoring.tsv` and
`data/metadata/enrichment/2026-10-04-physalia-surface/`. The preceding enriched
manifest remains preserved in `data/metadata/enrichment/2026-10-04/manifest.csv`.

Future depth summaries should show the source basis and retain collection
intervals. Sensitivity analyses can omit assigned depths; any comparison
between Physalia and deeper-collected taxa combines habitat with host identity
and study design. No depth-association test is established by this scoring step.

## Publication interval midpoints, 2026-10-12

The user subsequently requested a point value equal to the average of minimum
and maximum when no point depth is recorded. The publication metadata export
applies this to 14 libraries and labels the basis `interval_midpoint`, retaining
both bounds and the original depth text. Recorded point depths, including zero,
take precedence. NA19, WS5 and WS6 remain unresolved. This convention is applied
only in the [publication metadata](sample_metadata.md); the frozen manifest and
previously generated analyses are unchanged.

## Author-confirmed depths and event mean, 2026-10-21

Casey Dunn resolved NA19 to 296 m and WS5 (D325-SS12) to 246 m. For WS6
(D326-SS1), he supplied six animal collection depths and instructed use of
their arithmetic mean: 1971.53797 / 6 = 328.589661666... m, stored as
328.589662 m. The basis is `event_mean`, and bounds 295.8355–350.06003 m
retain the observed event range. The six source records are preserved in
[2026-10-21_ws6_collection_depths.tsv](../data/metadata/manual-review/2026-10-21_ws6_collection_depths.tsv).
This value is a representative event mean, not an interval midpoint or a measured
depth for the individual sequenced specimen. The collection-context reader
distinguishes it as `range_with_derived_mean` and uses the full range to assign
depth strata. This resolved the last numeric-depth gaps in the original 205
libraries. With the eight Bargmannia libraries, all 213 now have numeric depths.
Existing analyses and source snapshots retain their previous provenance;
the metadata export alone does not regenerate those results. The explicitly
authorized [expanded reporting refresh](expanded_reporting.md) incorporates
these reviewed values into current descriptive tables.
