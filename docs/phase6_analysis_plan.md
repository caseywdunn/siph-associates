# Phase-6 pre-specified analysis plan

Approved by Casey Dunn on 2026-09-29. It was fixed before any presence result
was summarized by host, region, study, or batch, as EXECUTION_PLAN.md Phase 6
requires ("predefine primary and sensitivity analyses before inspecting final
plots"). Settings are frozen in `config/phase6_analysis.json`. Every reported
table and figure must regenerate from the frozen long-form Phase-5 results
with the Phase-6 workflow. Any deviation from this plan is reported as a
deviation, with its reason.

## Design facts known at planning time (library metadata only)

- 205 libraries from 37 host species, 26 of them singletons. The replicated
  strata are *Physalia* (151 libraries, all Church 2025) and *Nanomia*
  (21 libraries, Ahuja 2026).
- 54 flowcells: 32 carry more than one host species, but only 3 carry more
  than one study, so study and batch are largely confounded.
- No blank or negative extraction controls exist.

## A. Evidence grades (per genome or vOTU, per library)

| Grade | Rule |
|---|---|
| Nominated | Phase-2 screen hit only |
| Trace | Bacteria 1–10% breadth; viruses 10–75% breadth |
| Validated | Phase-5 rule: bacteria ≥10% breadth and ≥100 reads; viruses ≥75% breadth |
| High confidence | Validated, plus same-library assembly support: a MAG of that species assembled from this library, or ≥10 kb of this library's Phase-3 contigs aligning to the genome at ≥95% identity (minimap2 `-x asm5`) |
| Probable contaminant (genome-level flag) | Validated presences concentrate on flowcells beyond host species and ocean region (test below) |

- **Reporting by host route.** Grades are reported by host route.
  Reference-free libraries were assembled from a 25 M-pair subsample, so they
  reach *high confidence* less often for technical reasons.
- **Contamination test.** For each genome with at least 3 validated presences:
  - The statistic is flowcell concentration, the sum over flowcells of
    (presences on flowcell)² ÷ total presences.
  - It is compared with 10,000 permutations of flowcell labels within host
    species × ocean-region strata (seed 20260929).
  - A genome is flagged at Benjamini–Hochberg FDR <0.05.
  - Strata with fewer than 2 flowcells do not contribute. A genome with no
    permutable stratum is reported as *untestable*, not clean.
  - Stratifying by region keeps real geographic associates, whose specimens
    were probably sequenced together, from being flagged.
- **Identity annotation.** Common reagent and skin genera (*Cutibacterium*,
  *Lawsonella*, *Bradyrhizobium*, …; listed in the config) are annotated
  separately. Contamination is never graded from identity alone.
- **Stated limitation.** Contamination spread evenly across batches (for
  example from handling) cannot be detected without blanks.

## B. Primary analyses

1. **Incidence and evidence-grade summaries** by host species and study, for
   every catalog genome and vOTU.
2. ***Physalia*** (151 libraries):
   - For each taxon validated in at least 10 libraries, fit a logistic mixed
     model: presence ~ log10(input pairs) + ocean region + (1 | flowcell)
     (lme4).
   - Region is tested by likelihood ratio against the model without region,
     with BH FDR across taxa.
   - Models are fitted only if region is identifiable apart from flowcell: at
     least 2 regions each on at least 2 flowcells, with flowcells spanning
     regions. Otherwise region patterns are reported descriptively, with the
     confounding stated.
   - Community level: Jaccard PERMANOVA on validated presence (vegan
     `adonis2`), ~ log10(input pairs) + ocean region, 9,999 permutations
     restricted within flowcell.
3. ***Nanomia*** (21 libraries): descriptive incidence by species, plus the
   same flowcell permutation test. No models.
4. **Broad host skim:** phylogeny-ordered descriptive summaries only, stating
   singleton species and study confounding explicitly.
5. **Cross-domain links:** phage–host pairs are reported only where supported
   by CRISPR spacers (catalog v1: 6 links). Other co-occurrence is described
   as association, never as infection.

## C. Sensitivity analyses

1. **Presence threshold.** Repeat B1–B3 with the bacterial presence threshold
   at 5% and at 20% breadth.
2. **Host handling.** Map 5 *P. physalis* and 5 *N. septata* libraries with
   full trimmed reads instead of host-depleted reads. They are a seeded draw
   (seed 20260930) across input-depth quintiles among assembly-eligible
   libraries. Compare presence calls per library.
3. **Capping.** Map 8 of the libraries above the 200 M-pair cap without
   capping. They are a seeded draw (seed 20261001) spanning host routes.
   Compare detections, especially trace-to-validated changes.
4. **Leave one study out** (no remapping). Count presences whose catalog
   genome was nominated only by one study's libraries. This does not capture
   read reallocation among the remaining genomes; a full rebuild and remap is
   run only if this shows heavy single-study dependence.

C2 and C3 need about 18 additional library mappings.

## Software

The statistical models and PERMANOVA run in a pinned R environment
(`envs/stats.yaml`: R with lme4 and vegan). Tabulation and permutation tests
run in the workflow's standard Python.

## Deviations

### 2026-09-30: per-taxon *Physalia* region test (approved by Casey Dunn)

- **What happened.** All 12 pre-specified logistic mixed models (B2) failed to
  converge: degenerate Hessians and singular fits. The cause is separation: a
  taxon is present in all or none of some region's libraries, so region
  effects are not estimable.
- **What is still reported.** The mixed models remain in the outputs
  (`physalia_models.tsv`), each labelled with its convergence status. None of
  their p-values is presented as a result.
- **Replacement test.** The reported per-taxon region test is now a
  permutation test that mirrors the pre-specified PERMANOVA design:
  - the statistic is the deviance gain from adding ocean region to a
    logistic model of presence on log10 input pairs;
  - it is compared with 9,999 permutations of region labels within flowcell
    (seed 20260930), with BH FDR across taxa;
  - output: `physalia_region_permutation.tsv`.
- **Why this test.** It asks the same question (region beyond depth, with
  batch accounted for) and does not depend on model estimates converging.
- **Alternative considered.** Weakly informative priors (`blme::bglmer`)
  would keep the model form, but they give less standard p-values.

### 2026-10-03: physical-flowcell correction and complete breadth sensitivity

A manuscript-readiness review found that the R regional analyses and Python
contamination test split `sequencing_batches` on `/`, although the manifest
uses `instrument:run:flowcell:lane`. Consequently, the original v1/v2 statistical
outputs treated lanes, and combinations of lanes, as separate flowcells. Those
batch-dependent results are superseded by the corrected run; no read mapping,
presence grading, assembly, or catalog reconstruction is required.

**Decision (approved by Casey Dunn, 2026-10-03).** Retain every library in
all descriptive results. Restrict flowcell-based inference to libraries
sequenced on one physical flowcell. A library spanning several lanes of the
same flowcell remains eligible; a library spanning different flowcells does
not. Do not encode overlapping flowcell sets as independent batch categories.
This restriction applies to the regional mixed models, regional permutation
tests, PERMANOVA, and contamination tests only. It does not apply to incidence,
evidence grades, eukaryote detections, assemblies, or comparative genomics.

`workflow/scripts/phase6_metadata.py` writes all 205 primary libraries to
`phase6_analysis/v2/primary/library_flowcells.tsv`, with physical-flowcell
identifiers, counts, inferential eligibility, and exclusion reasons. There are
136 single-flowcell libraries and 69 multi-flowcell libraries: respectively
122/29 in Church2025 (*Physalia*), 11/10 in Ahuja2026, and 3/30 in
Ahuja2024. Across studies, the genus *Nanomia* has 22 libraries (12 eligible
and 10 multi-flowcell), including one Ahuja2024 library in addition to the
21 Ahuja2026 libraries. The inferential
restriction changes the population to which the results apply; conclusions
about all libraries must therefore use the descriptive results.

For contamination tests, the minimum of three presences is assessed among
eligible libraries. Only host-species × region strata containing at least two
physical flowcells contribute to the concentration statistic and its
permutations, as specified in A. Results preserve total, eligible, permutable,
and excluded presence counts. An insufficient eligible count, absence of
permutable presences, or constant permutation statistic is explicitly
untestable; none is evidence of a clean sample or genome. The all-cohort and
*Nanomia*-only tests are reported separately, with BH adjustment within each
test family and threshold. Sorted strata make the seeded randomization
independent of Python hash/set iteration order.

The previous breadth sensitivity table contained only aggregate totals,
although C1 required repeating B1–B3. The corrected workflow repeats incidence,
evidence summaries, contamination tests (including the *Nanomia* subset),
*Physalia* mixed models, per-taxon regional permutation tests, and community
PERMANOVA at 5% and 20% bacterial breadth. Results live under
`sensitivity/breadth_5pct/` and `sensitivity/breadth_20pct/`; primary 10% outputs
retain their original paths. Each threshold selects taxa with at least ten
presences among the eligible *Physalia* libraries, so tested families may
change. All per-taxon tests use the same 9,999 seeded permutations; computational
caching and parallel execution do not alter the statistic.

PERMANOVA uses eligible libraries with at least one bacterial/archaeal
presence. Empty communities are excluded because binary Jaccard distance is
undefined for two empty communities; each threshold reports the eligible,
included, and empty-community counts. Regional identifiability is checked
again on this nonempty subset. Community and per-taxon regional results remain
observational associations, with host species, region, study, and technical
sampling limitations retained in their interpretation.
