# Host sampling and descriptive detection summaries

This reporting step joins reviewed specimen metadata, accepted detections and
selected assembly membership. It supplies the host denominators and descriptive
counts used in the paper. Current input paths are selected by
[reporting_sources.json](../config/reporting_sources.json); no sequence analysis,
reference change or hypothesis test is run by this script.

The cohort contains 151 *Physalia*, 22 *Nanomia* and 40 other siphonophores.
Among the ten *N. bijuga*, four have accepted bacterial detections. There are
207 specimens assigned to 37 named taxa, with six other specimens retaining
broader or provisional identifications. Current identities follow the
[publication host-name table](../config/publication_host_names.tsv) and
[specimen metadata](../figures/cohort/sample_metadata.tsv).

The frozen analysis identifiers and names reported by the source studies are
retained. This includes the five *Physalia mikazuki* formerly labelled sp1/B2
and the [reviewed Nanomia identities](sample_metadata.md#nanomia-identifications).
Unresolved host labels are kept separately and excluded from named-taxon counts.

## Rebuild the summaries

First produce current cohort and collection-context tables through the
[reviewed reporting workflow](viral_reporting_review.md). From the analysis
repository root with Python 3.11 or later (standard library only), the bounded
summary can also be run independently:

```bash
python3 scripts/summarize_dalila_revision.py
```

Outputs in [figures/dalila_revision](../figures/dalila_revision/) include generated
manuscript values, all-host detection denominators, regional *Physalia* counts,
recorded life-stage entries and provenance with input/output and code checksums.
The `expanded_supporting_tables` and reviewed reporting stages also run this
producer. A missing life stage remains unknown; it is not interpreted as an
adult. These summaries do not establish independent age, depth, host or regional
effects.
