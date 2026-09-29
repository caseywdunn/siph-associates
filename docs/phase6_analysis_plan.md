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
