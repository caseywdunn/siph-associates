# Catalog interpretation for reporting

These manifests derive from the frozen v3 mapping catalog using
[`config/viral_catalog_classifications.tsv`](../../../config/viral_catalog_classifications.tsv).
They preserve all target identifiers and sequence checksums. Three
Negarnaviricota targets are reported as possible endogenous elements because
their evidence comes from DNA sequencing; integration has not been demonstrated.

The bacterial and MAG manifests are unchanged copies included for a coherent
reporting source selection. `provenance.json` records original input hashes and
identical reused-output hashes. `summary.json` records the resulting class
counts, and `reclassified_detections.tsv` lists the seven accepted observations
whose interpretation changed. Original mappings and the original catalog
manifests remain unchanged.

Rebuild through [`workflow/coauthor_reporting.smk`](../../../workflow/coauthor_reporting.smk),
as described in the [reporting guide](../../../docs/viral_reporting_review.md).
