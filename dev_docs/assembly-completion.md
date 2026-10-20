# Assembly completion and mapping preparation, 2026-10-20

User removed the retained-read-count cutoff: attempt every library, accept
no-contig outcomes, prepare the expanded cohort, and do not run full-cohort
mapping. Additional MAG/viral identity and contamination QC precedes catalog
approval and any mapping authorization.

The current configuration is `config/assembly_completion.json`; its generated
manifests select 12 new assembly attempts and exactly 213 cohort libraries.
Original 175k cutoff configurations/membership remain historical records.
The selected cohort reuses 193 original assemblies, replaces the two original
Halistemma full-read assemblies, includes the ten previously excluded inputs,
and adds eight completed Bargmannia assemblies under their new species names.

Execution **13631472** completed all 12 assemblies, marker searches, viral
searches and CheckV. Binning failed while archiving because the new entry point
supplied a relative scratch path to a shared rule that changes directory before
creating the archive. No scientific tool parameters were changed to fix this:
the new entry point now supplies an absolute scratch path. Source snapshots and
failed logs remain in `data/results/assembly_completion/2026-10-20/executions/`.

The recovery dry run scheduled only 12 binning steps, 12 CheckM2 steps, two
summaries/checks and an aggregate target. **13632101** completed that scope at
17:28 EDT, exit 0; completed assemblies and marker/viral outputs were reused.
The stage validator reports PASS, 12 assembled, zero excluded, and zero errors.

The 12 assemblies all produced contigs. Early inspected examples: CWD12 has
322 contigs / 470,091 bases, and YPM-IZ-111013 (88,861 input pairs) has 18 contigs
/ 32,844 bases. These are assembly yields, not validated associate assignments.

Total new assembly yield is 4,196 contigs / 8,607,764 bases, two raw bins and
38 geNomad viral candidates. One bin, from YPM-IZ-111016 (145,748 retained
pairs), passes the existing CheckM2 screening criteria. Its domain, taxonomic
identity and contamination remain to be reviewed. Thus the old input cutoff
would have excluded a genome candidate; assembly yield was not determined by
that cutoff.

Validation includes unique replacement selection, inclusion of tiny/zero-count
inputs in the plan, acceptance of empty assembly outcomes by the shared validator,
preservation of explicit mapping holds, candidate sequence/ID checks, and
prevention of substituting 25M-pair assembly subsets for full mapping inputs.
Five focused tests and Ruff passed. All 12 selected input records matched the
completed mapping/trim provenance and retained FASTQ sizes.

The mapping-input inventory has 163 host-depleted libraries and 50 full capped
libraries. Files are retained for 171; 42 original reference-free libraries
require trimmed-read regeneration when mapping is authorized. Preparation does
not regenerate them. The hold now requires explicit authorization even after
all assembly and preparation jobs finish.

Preparation **13632477** completed all four steps at 17:31 EDT. It verified
213 selected assembly checksums and staged 1,917 links to completed products.
The review bundle contains 131 CheckM2-screened bins (127 reused, one from
assembly completion, three from Bargmannia), including the four sequence-matched
known eukaryotic controls. There are 34,566 candidate contigs. All gene-domain
and contig-taxonomy reviews remain pending; no catalog is approved.

Final dry runs report both workflows up to date. Both generated rule graphs
match their workflows. Mapping remains unauthorized and was not launched.
The candidate table and preparation summary are in `figures/catalog_qc/`;
the source and mapping-input inventories are under
`data/results/post_expansion/2026-10-20/`.

Canonical commands and outputs: [completion/preparation guide](../docs/assembly_completion.md).
Required deeper QC: [catalog QC plan](../docs/catalog_qc.md).
