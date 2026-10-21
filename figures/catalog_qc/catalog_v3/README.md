# Curated catalog v3

These tables record the expanded catalogs after MAG domain review, foreign
fragment removal, updated CheckM2 estimates and viral boundary curation.
[Methods and execution records](../../../docs/catalog_qc.md#curated-mapping-targets)
describe the decisions and bounded Snakemake workflows.

| Catalog component | Count |
| --- | ---: |
| Retained bacterial MAGs | 104 |
| MAG clusters | 32 |
| Pooled microbial targets | 478 |
| MAG representatives in the pooled targets | 24 |
| Public reference representatives in the pooled targets | 454 |
| Additional fixed decoys | 20 |
| Associate-virus clusters | 407 |
| Endogenous-virus candidate clusters | 157 |
| CRISPR host links | 6 |

Eight MAG clusters have public reference representatives in the pooled catalog.
MAGSP0033 is the new Tyloplasma cluster from YPM-IZ-111016; its pooled target is
BAC00459. The false archaeal identifiers MAGSP0026 and BAC00006 remain retired.

- [Readiness report](catalog_ready.json): validation result and mapping hold at
  the original pre-mapping checkpoint; subsequent authorization is recorded separately.
- [File inventory](inventory.tsv): paths, sizes and SHA-256 checksums for the
  final FASTAs, indexes and review tables.
- [Identifier crosswalk](catalog_id_crosswalk.tsv): old/new cluster membership
  relationships, retained identifiers, additions and exclusions.
- [MAG catalog](mag_catalog.tsv): all accepted MAGs and their clusters.
- [Microbial manifest](bacterial_catalog.manifest.tsv): pooled targets,
  reference exclusions, clustering and decoys.
- [Associate-virus manifest](viral_associate.manifest.tsv) and
  [endogenous-virus manifest](viral_endogenous_candidate.manifest.tsv): viral
  representatives and sequence checksums.
- [CRISPR host links](host_links.tsv): sequence-based links between catalog
  genomes and viruses.

The FASTAs and indexes remain on the cluster at
`data/results/post_expansion/2026-10-20/curated/phase4_catalog/v3/`.
The microbial and combined viral FASTAs support separate competitive alignments
in the same cohort-remapping round. This record was created before mapping;
the subsequent [authorization and run guide](../../../docs/cohort_remapping.md)
records the completed 213-library remapping. The
[current reporting guide](../../../docs/expanded_reporting.md) identifies v3
grades, statistics and publication exports. This immutable catalog receipt
retains the hold state that applied when it was created.
