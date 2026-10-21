# Catalog correction and current statistics

`mag_catalog.tsv`, `bacterial_catalog.manifest.tsv`, the gene/marker domain
tables, `summary.json`, `values.tex` and `provenance.json` preserve the original
false-archaeal correction. They describe that historical v2 review, not the
current expanded genome inventory.

The current inventory is in [catalog_v3](../catalog_qc/catalog_v3/):
[MAGs](../catalog_qc/catalog_v3/mag_catalog.tsv),
[microbial representatives](../catalog_qc/catalog_v3/bacterial_catalog.manifest.tsv)
and the [catalog crosswalk](../catalog_qc/catalog_v3/catalog_id_crosswalk.tsv).
The incorrectly assigned `MAGSP0026` remains excluded.

The `statistics.*` and `threshold_*.tex` files in this directory are refreshed
from the current v3 grades and geographic/batch tests. Their provenance file
records the specific inputs. They are produced by `expanded_supporting_tables`
in the [reporting workflow](../../docs/expanded_reporting.md).
