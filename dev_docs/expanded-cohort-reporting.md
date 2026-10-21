# Expanded-cohort reporting, 2026-10-21

Scope authorized by Casey Dunn: review the completed v3 mappings, regenerate
downstream associate summaries, and bring analysis/workflow documentation up to
date before editing the manuscript. Source checkpoint: `44ade17`.

The mapping controller `13646027` completed at 14:12 EDT on October 21 after
939/939 steps. Its final record reports PASS for 213 libraries, 498 microbial
targets including 20 decoys, 564 viral targets, and all 42 regenerated-read
checks. Execution records remain under
`data/results/post_expansion/2026-10-20/remapping/executions/13646027/`.

## Work and validation

- Preserve the current publication outputs and their source records as the
  comparison baseline. Reuse the completed mappings, assemblies and screens.
- Use the reviewed v3 catalog and selected assemblies for same-library assembly
  support. Retain the locked presence thresholds and source-MAG/vOTU support.
- Recompute grades, incidence, physical-flowcell tests, geographic tests and
  threshold sensitivity on the expanded cohort. Distinguish any historical
  control analyses from newly computed v3 results; do not silently mix catalogs.
- Incorporate the eight new libraries into publication metadata and SSU evidence,
  including appropriate host-reference competition and source verification.
- Audit catalog-dependent phylogeny/function inputs and retain unaffected
  computations only when their sequence identities and evidence still apply.
- Rebuild analysis-repository tables and figures, compare old/new results using
  stable sample and catalog identities, and record biological changes.
- Update active README, workflow guides, graphs and execution records; preserve
  dated historical audits. Check links, graph freshness, bounded dry runs and
  relevant tests. Manuscript edits are outside this task.

Status: complete. Final allocation 13687331 passed 14/14 reporting steps on
2026-10-21 at 16:26 EDT. Both validation sentinels report PASS with zero errors.
The completed publication dry run reports nothing to do. The chronological
execution notes below retain intermediate states and recovery decisions.

Execution started with controllers 13685942 (statistics) and 13686138 (local
SSU evidence). The latter is being consolidated into one Slurm allocation to
avoid the 200-job/hour submission limit; scientific inputs and rules are unchanged.
Automatic approval review rejected a target that could send new assembled SSUs
to NCBI BLAST. Only local evidence work has been launched; sequence-identical
cached verification remains reusable without a new external query.

Local continuation 13686930 completed 229/229 steps, including all 213 assembled
SSU classifications, expanded screen evidence and Phase-6 validation (PASS,
zero errors). Normalization accepts both physical-flowcell IDs and lane-qualified
IDs; the current inference pool is 144 libraries, including 122 Physalia.
All 65 assembled non-host queries match the old sequence cache exactly.
Verification is now explicitly cache-only, failing before network access for any
novel sequence; no export approval is needed. Array 13686929 continues six
read-classification batches; the final unbatched target completes the seventh
batch and aggregation. Pending allocation 13686507 was cancelled before starting.

Completed checks before the final exports:

- All 56 focused tests passed (metadata, collection context, phase-6 statistics,
  cache-only verification, eukaryote reporting, host-depletion reporting and
  fish classification).
- A fresh `workflow/cohort_remapping.smk remapping --dry-run --rerun-triggers
  mtime` reports nothing to do. The expanded publication dry run schedules
  only unfinished SSU classifications and downstream reports.
- Fish reports were merged across all 213 libraries: 9,577,330,361 pairs,
  32 strong pairs in eight Physalia libraries, no new strong Bargmannia hits.
- v3 has 568 microbial detections: all 557 original detections remain, with
  nine in new libraries and two Tyloplasma detections in original libraries.
  Four original detection grades change with updated assembly support.
- Associate-virus detections total 631. After catalog-ID translation, 13 old
  detections are lost on six non-retained targets, 15 are gained in original
  libraries, and 99 occur in the eight added libraries. Of these additions,
  86 occur in formalin-preserved Bargmannia elongata T676-SS3.
- Public documentation now distinguishes the current v3 route from dated
  original design/calibration records. The seven unchanged focal genome
  analyses retain their scope; the new Tyloplasma genome is not represented in
  their tree or functional comparison.

Array 13686929 task 2 started after the two retained long-running read jobs
completed. To avoid serial waiting, allocation 13688086 processes the final
33-library subset (27 unfinished tables when submitted), with a dependency on
completed array task 6. Its exact targets are recorded in
`executions/final_ssu_batch.json`; they are the final 33 sorted aggregation
inputs, disjoint from the six 30-input array batches. Final publication job
13687331 now waits on both the array and 13688086. At most six 32-core SSU
allocations run together. Readiness still requires final exports, biological
comparison, figure inspection and the remaining publication graph check.

The v3 microbial total includes 567 bacterial detections and one archaeal
reference detection. D666-D10 maps to Nitrosopelagicus reference
GCA_013390375.1/BAC00450 (breadth 0.131901, 1,083 retained reads, zero passing
assembly-support bases); it is trace at the 20% breadth sensitivity cutoff.
This is not an archaeal MAG. The reporting guard that had assumed all archaeal
detections were MAG-only was replaced with explicit MAG/reference counts.
The previously rejected MAGSP0026 remains excluded. Domain-specific cohort
summary counts now distinguish bacterial, archaeal and combined microbial totals.

The two MAGSP0006/BAC00362 calls downgraded from high confidence to validated
are YPM-IZ-110467 and YPM-IZ-110468. The representative lost four host/eukaryotic
contigs (12,570 bp) during the accepted QC curation; passing assembly support
decreased from 14,380/15,825 to 2,698/4,029 bp. Both retain read-supported
presence. The contig decisions are in `figures/catalog_qc/contig_decisions.tsv`.

All 35 graph scopes other than the bounded publication graph passed the
freshness check, and all 36 SVG files parsed. The expanded publication graph
was visually inspected. All 46 then-changed Markdown files had valid local
file links. The bounded publication graph awaits the final eukaryote table.

After four array tasks finished, batch 2 still held most of the unfinished
large SSU inputs. Its Snakemake 7.24.0 scheduler source was checked: SIGTERM to
the controller invokes `exit_gracefully`, which stops new scheduling and waits
for active jobs. Only controller PID 2476123 was signalled in allocation
13686929_2, and its log confirmed this behavior. Four active searches were
left running; all completed files were retained. The 14 unstarted samples were
partitioned by compressed input size into three disjoint local allocations:
13688502 (4 libraries), 13688503 (4), and 13688504 (6). The exact membership,
active/completed exclusions and final dependencies are recorded in
`executions/ssu_batch2_redistribution.json`.

The publication dependency is now `afterany:13686929` plus successful completion
of 13688086/13688502/13688503/13688504. The intentionally drained batch exits
nonzero, so `afterany` is deliberate; final Snakemake dependencies and validators
still require every library output. This rescheduling preserves at most six
32-core SSU allocations, changes no search command and repeats no finished or
active search. Source archives remain attached to each execution.

The completed T676-SS3 read table exposes a specimen-specific background issue:
288 insect and 71 plant pairs at ≥97% identity exceed the locked 37-pair
threshold; the previous cohort maxima remain 36 and 13. A generated flag table
and descriptive sensitivity summary omit flagged libraries across each domain
without altering primary grades or retuning the threshold. The predefined
background classification is imported from the shared SSU reporting module.
A focused test checks identity-band aggregation, the 36/37 boundary and exclusion
of host/human assignments. This SSU result does not prove the 86 viral calls in
the same specimen false.

The supplied `tmp/zlag150.pdf`, molecular-methods section, identifies T676-SS3
as a formalin-preserved nectophore extracted with a modified Qiagen protocol
(washing, autoclaving and extended digestion); the other seven specimens used
frozen tissue and an E.Z.N.A. kit. It describes common library preparation at
YCGA, without reporting whole-genome amplification or identifying a cause for
the background signal. No new mechanistic attribution is made. Publication
allocation 13687331 was briefly held before starting to include this review,
then released with the final SSU dependencies intact.

## Completed comparison and handoff

All 213 competitive read tables completed; the last allocation 13688502 ended
successfully at 16:25:56 EDT. Final reporting 13687331 completed at 16:26:24,
exit 0. Its source archive, execution log and finish record are under
`executions/13687331/`. Phase-6 and eukaryote validators both report zero errors.
The eukaryote result has 67 validated evidence rows, 64 nonredundant detections,
40 positive libraries, 17 groups and 42 high-confidence detections. All original
61 detections and their grades remain; three read-only groups are added in
T676-SS3. All 65 assembled non-host sequences used cached verification, with
195 retained hit rows and no new query.

The generated comparison confirms zero changes to the original 205 libraries'
reviewed metadata. All original microbial presences remain. Excluding the only
background-flagged library gives 567 microbial calls (139 libraries), 545 viral
calls (126 libraries), and 61 eukaryotic calls (39 libraries). Its 86 viral
targets are unique to that specimen at the accepted threshold. These results
qualify interpretation without changing primary evidence grades.

Final checks: 57 focused tests passed; new helpers passed Ruff formatting and
lint; shell launchers passed syntax checks; Markdown links resolved. All 36
graph scopes have current DOT/SVG checks, including the final bounded and
expanded publication scopes. The refreshed microbial/eukaryotic overview and
three collection-context figures were visually inspected. The unused overview
concept sketch retains explicit photograph placeholders; the manuscript's
Figure 1 is maintained separately. The remapping and final publication dry
runs both report no work. No costly clean-state sequence rerun was performed.

The [comparison](../figures/expanded_review/README.md) and
[manuscript handoff](expanded-manuscript-handoff.md) are current. No manuscript
files were changed. Remaining biological uncertainty includes T676-SS3's
background, the read-only archaeal reference match and previously documented
metadata source conflicts; successful execution does not resolve those issues.

A final reporting-only refresh, allocation 13690006, passed 5/5 steps at
16:29:16 EDT and corrected two descriptive
strings in the collection-context provenance: interval strata use the full
bounds while derived point depths remain recorded, and comparisons are labelled
descriptive. It repeats no sequence search or statistical inference. Numeric
results are unchanged.
