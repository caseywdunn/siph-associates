# Removed gene-content analysis, 2026-10-24

The author removed predicted nutrient use, metabolism, defense genes and nutritional co-occurrence from the paper and maintained analysis code. Historical source is retained in Git; untracked computational outputs remain.

The focused phylogeny now reads its study/external labels from the retained tree-label table and GTDB r220 taxonomy from `config/mycoplasmatales_reference_taxonomy.tsv`. The latter retains only the genome identifiers and classifications needed for plotting from the previous summary, with no functional annotations.

- Source: `data/results/mycoplasmatales_gene_content/summary/genome_stats.tsv`
- Source SHA-256: `87c8a246b3e8123ceedfd33646a82dfeb819e2e2157de446ea0362bc197dd189`
- Base revision: `2a4c03da69e8e432cecabca4e637f740261a9bcd`
- Taxonomy rows: 142
- MAGSP0031 is labeled a sister lineage of clade B. Tree inference is unchanged.
- Execution and workflow checks are pending.
