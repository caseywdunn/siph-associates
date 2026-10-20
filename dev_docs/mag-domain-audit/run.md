# MAGSP0026 domain review, 2026-10-18

## Completed review and correction, 2026-10-19

Audit job 13526694 completed all 15 steps with exit code 0; all 15 reference
controls passed. nr assigned 2,518 of 3,615 representative proteins to Eukaryota
and none unambiguously to Archaea. All 29 representative archaeal-marker hits
match Perkinsela proteins. The other three bins show the same pattern.
These results supersede the initial pending assessment below.

`config/catalog_domain_exclusions.tsv` excludes BAC00006/MAGSP0026 from
prokaryotic reporting. Original catalog and sequence files remain archived;
reviewed catalog tables are in `figures/catalog_review/`. The other bacterial
grade rows are exactly unchanged. Removal of 205 matrix rows (one target in
each library) removes 14 previously accepted detections and four retained MAGs.
The retained prokaryotic set is 103 bacterial MAGs in 31 clusters.

Correction job 13567655 completed successfully. Its dry run and preserved source
are under `data/results/catalog_domain_review/2026-10-19/`; the pre-correction
report tables are preserved in `pre_correction.tar.gz`. It reran only grading,
descriptive exports and statistical tests, with no sequence mapping or searches.
The Phase-6 validator and targeted removal checks passed. Supplementary reporting
then regenerated the unchanged focal bacterial summaries and the corrected
overview figure. `scripts/export_domain_review_statistics.py` generates affected
manuscript values and threshold-table rows. Manuscript correction uses these
exports and retains the bins for eukaryotic follow-up without automatically
adding genome-mapping evidence to the separately calibrated SSU detection counts.

Scope: MAGSP0026 and its three related bins. Requested after an archaeal
specialist reported protist/trypanosome-like genes and no archaeal genes.

The FASTA provided to the specialist matches the original archived bin,
candidate and v1/current representatives, SHA-256
`2177988a884aeea1cda368e31d2586c3956b349e69172746472db71ade18f712`.
The original GTDB-Tk r220 summary reports RED-based Thermoproteota placement,
40.11% marker-alignment occupancy and no reported reference-genome ANI match.
Only non-diagnostic 5S rRNA predictions occur in this bin. Initial review did
not establish its true domain; the specialist's gene searches have not yet
been independently repeated here.

Worktree: `data/results/worktrees/mag-domain-audit`, branch `mag-domain-audit`,
based on `master` at `02a1c83`. Existing fish work and its source files remain
in the original worktree. No sequence-processing rule from the original
analysis is included in the new workflow.

Methods and commands: [domain audit](../../docs/mag_domain_audit.md).
Validation: five focused tests pass; Python lint and formatting pass; shell
syntax passes. The dry run contains 15 jobs, all within the new audit scope.
Production searches are gated on known-reference sequence/taxonomy controls.
Release/download checks found Swiss-Prot 2026_03 and an accessible NCBI taxonomy
archive. An initially incorrect release-notes URL was corrected before launch.

Submitted with `sbatch --parsable scripts/mag_domain_audit/batch.sh`: job
**13526694**, confirmed RUNNING on `r814u03n03` with 16 CPUs and 64 GB.
It executes the preserved source snapshot at
`data/results/mag_domain_audit/2026-10-18/executions/13526694/source/` in the
original repository. The source archive, source checksum, environment records,
execution log, and eventual exit code are retained beside it. Both tool-version
checks succeeded (DIAMOND 2.2.0 and GTDB-Tk 2.4.1).

Rule-graph generation and its freshness check passed. Production controls and
scientific assessment remain pending while the job runs. Fish jobs were not
cancelled or resubmitted.

Initial execution check: Swiss-Prot 2026_03 and NCBI taxonomy downloaded;
Swiss-Prot index contains 575,748 proteins. GTDB marker identification and gene
collection completed: 12,463 predicted proteins across the four bins, including
3,615 in the MAGSP0026 representative. The marker-hit table retains 26, 29, 19
and 26 hits respectively in the manifest order, including multiple-copy markers;
these counts do not establish archaeal identity. Reference-protein controls were
running at the initial status check.
