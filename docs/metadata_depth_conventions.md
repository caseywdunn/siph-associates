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
