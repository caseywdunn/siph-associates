# Viral interpretation and host-name reporting

The selected reporting tables apply reviewed biological labels to completed
mapping evidence. Raw mappings, mapping target identifiers, reference sequences,
and evidence thresholds remain those of catalog v3. The current source paths
are selected in [reporting_sources.json](../config/reporting_sources.json).

[viral_catalog_classifications.tsv](../config/viral_catalog_classifications.tsv)
classifies VOTU00144, VOTU00165 and VOTU00242 as endogenous candidates. These
Negarnaviricota assignments came from DNA libraries; they are possible endogenous
elements, and their integration into a genome has not been demonstrated.
VOTU00242 has seven accepted detections, three with assembly support. The other
two targets have no accepted detections. The target IDs are retained so every
observation remains traceable to the original competitive mapping.

The derived [catalog manifests](../figures/coauthor_review/catalog) contain
404 associate groups and 160 endogenous candidates. Of the associate groups,
391 are detected in 122 specimens, giving 624 detections and 439 with assembly
support. The original [mapping catalog](../figures/catalog_qc/catalog_v3)
remains the record of the references used for alignment.

Host labels come from the publication metadata overlay. The 22 *Nanomia*
comprise ten *N. bijuga*, ten *N. septata*, one *N. cara* (WS4) and one *Nanomia*
sp. 1 (NA29). [Voucher-level evidence](sample_metadata.md#nanomia-identifications)
resolves the earlier pooled labels; Namrata Ahuja separately confirmed NA37 as
*N. bijuga*. These labels do not change library IDs or processing routes.
Inconsistent life-stage labels remain in the source metadata and are excluded
from biological summaries and tests.

## Rebuild from completed evidence

Use the existing Snakemake environment from the repository root:

```bash
export XDG_CACHE_HOME=/tmp/siph-snakemake-cache
snakemake -s workflow/coauthor_reporting.smk --cores 2 --dry-run
snakemake -s workflow/coauthor_reporting.smk --cores 2
```

This entry point contains only metadata, interpretation, statistical summary and
figure rules. Missing completed input tables cause a failure; it has no FASTQ,
assembly or read-mapping producers. Flowcell concentration tests and incidence
summaries are recomputed at all three bacterial breadth thresholds because
reviewed *Nanomia* identifications change host strata. Regional tests concern *Physalia* and reuse the
validated completed outputs after checking that its metadata and bacterial
evidence are unchanged. Input and reused-output hashes are recorded in
`data/results/coauthor_review/2026-10-22/prepared.json`.

![Reporting workflow](workflow/graphs/coauthor_reporting.svg)

Validation checks the complete sample-by-target evidence matrices, all three
breadth thresholds, physical-flowcell eligibility, and unchanged mapping
measurements and grades. The selected cohort, collection-context, overview,
annotation/co-occurrence, background-control, screening and fish summaries are
then regenerated. No sequence analyses are run.
