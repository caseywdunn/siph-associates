# Museum metadata comparison, October 17, 2026

Reviewed 148 exact YPM voucher matches in the fresh
[LUX/GBIF audit](../2026-10-17-public-vouchers/summary.json) and 106 public locality
records linked from those specimens. The remaining 20 vouchers were still not
discoverable in LUX or GBIF. The locality links were initially taken from the
October 12 audit and verified against every October 17 object record before use.
The newly retrieved locality records and request receipts are retained in
`data/sources/metadata_enrichment/2026-10-17-lux/`.

`museum_metadata.tsv` retains the museum labels, field/personal identifiers,
locality descriptions, depth text, and mapped locality geometry for each voucher.
`field_comparisons.tsv` records agreements, differences, and comparison limits.
`review_needed.tsv` lists 14 date disagreements, 10 coordinate discrepancies,
2 depth disagreements, 6 alternate field labels, and 1 inconsistent country.
These categories can overlap within a specimen. Different identifiers with
different roles are retained in the full comparison rather than all being called
errors. An agreement does not establish that two sources are independent.

All 443 rows of `proposed_updates.tsv` were accepted and appended to the curation
overlay. Of these, 221 add museum corroboration for unchanged identifiers; 222
change publication metadata cells: 49 site descriptions, 49 location refinements,
26 country/territory entries, 22 collection methods/platforms, 38 depth-bound
cells for 19 specimens, and 38 human-facing notes. No existing overlay rows needed
replacement. The prior table and overlay are retained as `sample_metadata_before.tsv`
and `overlay_before.tsv`; `accepted_changes.json` records the accepted changes.

Existing point coordinates, point depths, collection dates, tissue descriptions,
and identifiers are retained. Conflicting museum values are explained in notes.
The CWD131 country is withheld because its museum locality combines Gulf of
California with USA. New site/country details are withheld where the study point
falls outside the museum locality envelope. A location polygon is context, not a
measured specimen coordinate; no polygon centres are imputed. A method-only
locality label does not establish a geographic site. For Johnson Sea Link II,
the explicit platform name is retained without adopting the catalog's vehicle
classification. Museum field/personal numbers are not used to fill missing
tissue labels. Museum internal locality numbers are not field collecting-event
identifiers.

Depths are parsed only when units are explicit. Compatible museum bounds augment
an existing study point depth; they do not replace it with a midpoint. The study's
0–20 m interval for CWD16 and 1299 m point for NA8 are retained despite museum
values of 0–50 m and 1289.68 m. The public NA19 locality does not report depth,
so it cannot resolve the prior 296/297 m source disagreement. No previously
missing required fields were resolved; the new information is collection context
and source verification.

Python 3.11 was used; the retrieval and review scripts use only the standard
library. Fresh retrieval uses a new dated output directory:

```bash
python3 scripts/check_public_ypm_vouchers.py \
  --output-dir data/metadata/enrichment/2026-10-17-public-vouchers
python3 scripts/fetch_lux_localities.py \
  --audit-dir data/metadata/enrichment/2026-10-17-public-vouchers \
  --output-dir data/sources/metadata_enrichment/2026-10-17-lux
```

Reproduce the retained comparison offline against the pre-curation snapshots:

```bash
python3 scripts/review_lux_metadata.py \
  --audit-dir data/metadata/enrichment/2026-10-17-public-vouchers \
  --localities-dir data/sources/metadata_enrichment/2026-10-17-lux \
  --metadata data/metadata/enrichment/2026-10-17-lux/sample_metadata_before.tsv \
  --overlay data/metadata/enrichment/2026-10-17-lux/overlay_before.tsv \
  --output-dir /tmp/lux-metadata-review
```

After reviewing a proposal and updating the durable overlay, the metadata-only
export is `python3 scripts/export_sample_metadata.py`. It does not refresh figures
or execute sequence analyses. Frozen input checksums and modification times are
recorded in `protected_inputs_before.json` and checked after export.
