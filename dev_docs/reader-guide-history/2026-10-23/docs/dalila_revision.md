# Dalila manuscript revision, 2026-10-11

The opening record below describes the October 11 revision. Its exclusion of
the Halistemma reference was superseded by the subsequently authorized
[Halistemma analyses](halistemma.md). The same descriptive reporting script
now consumes the 213-library metadata, v3 detections and selected assemblies
through `config/reporting_sources.json`; `expanded_supporting_tables` in the
[current reporting workflow](expanded_reporting.md) rebuilds these exports.

Reporting-only changes on source revision 00a629cd917440a97a8a5562b70e6837a5115d83.
No sequence analysis, reference update, catalog change or hypothesis-test rerun.
The user explicitly excluded incorporating the new Halistemma reference.

`config/publication_host_names.tsv` maps the five frozen Physalia sp1 / B2
specimens to Physalia mikazuki. Cluster correspondence: Yongstar et al. 2025,
10.3389/fmars.2025.1653958; available species description: Yongstar et al. 2026,
10.50974/0002008267. The original manifest and species_as_published are retained.
The overlay is applied once in figures/build_cohort.py; downstream display
summaries consume its sample_metadata.tsv. Saved detections and statistical
outputs are unchanged. The tree legend now identifies bacterial MAGs explicitly.

Run from the analysis root with the existing Snakemake Python environment:

```sh
python figures/build_cohort.py
python scripts/analyze_collection_context.py
python figures/plot_collection_context.py
python figures/plot_associate_overview.py
python workflow/scripts/plot_mycoplasmatales_figures.py
python scripts/summarize_mdb_revision.py
python scripts/summarize_dalila_revision.py
```

The final script joins reviewed sample metadata and collection-context detections,
producing generated manuscript values, all-host detection denominators, all-Physalia
regional denominators for clade A, and the recorded life-stage entries. Unresolved
host labels remain separate and are excluded from named-species denominators.
Missing life-stage entries are unknown, not adults. These are descriptive joins,
not new tests or evidence of independent age, depth, host or geographic effects.
Its provenance records input/output checksums and source code identity.
