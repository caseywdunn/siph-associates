# Handoff for the expanded-cohort manuscript update

Prepared 2026-10-21. Manuscript files have only been inspected, not edited.
**Status: analysis exports complete and validated; ready for manuscript
reconciliation with the T676-SS3 qualification below.** The current
execution state is in [expanded-cohort-reporting.md](expanded-cohort-reporting.md);
the reader-facing route is [expanded reporting](../docs/expanded_reporting.md).

## Selected scope and results

- Four studies, 213 libraries: Church2025 151, Ahuja2024 33, Ahuja2026 21,
  Manko2026 8. There are 151 Physalia, 22 Nanomia and 40 other libraries,
  with 42 current host labels. Bargmannia includes the published B. pyrosa
  and B. marlena names; the NCBI labels and collecting-event source conflicts
  remain separately documented.
- All 213 selected libraries yielded contigs after an assembly attempt.
  No retained-read minimum excludes a library. The selected inventory combines
  193 original assemblies, ten formerly excluded low-input assemblies, two
  Halistemma replacements and eight Bargmannia assemblies.
- Host depletion covers 151 Physalia, ten Nanomia septata and two Halistemma
  rubrum libraries. Current full-cohort mapping uses all retained depleted
  reads and up to 200 million raw pairs before quality filtering elsewhere.
- The curated inventory has 104 bacterial MAGs in 32 clusters, 55 near-complete.
  Catalog v3 contains 478 microbial representatives plus 20 decoys, and 407
  associate plus 157 endogenous-candidate viral representatives. Use the
  [v3 inventory](../figures/catalog_qc/catalog_v3/), not the historical
  false-archaeal-correction inventory in `figures/catalog_review/`.
- Primary microbial calls total 568 (567 bacterial and one archaeal), across
  140 libraries and 189 targets; 303 have assembly support. Every original
  microbial detection remains. Nine calls occur in Bargmannia and two new
  Tyloplasma calls occur in Physalia. Four original evidence grades change.
- The new Tyloplasma representative is MAGSP0033/BAC00459, recovered from
  YPM-IZ-111016 and also detected in YPM-IZ-111017. It belongs to
  Hepatoplasmataceae. The seven focal Metamycoplasmataceae genomes are
  sequence-identical to their previous versions, so their tree, annotations
  and external cnidarian comparisons are retained. Do not describe that
  seven-genome comparison as covering all newly recovered Mycoplasmatales.
- The eight Physalia source assemblies for the Waki marker comparison are
  byte-identical to their selected expanded-cohort assemblies. Retain that
  comparison and its eight-specimen scope; its 110931 host label is now
  Physalia mikazuki in descriptive documentation.
- The archaeal signal is a read-only Nitrosopelagicus reference match in
  Bargmannia stenotes D666-D10: BAC00450/GCA_013390375.1, 13.1901% breadth,
  1,083 retained reads, no passing assembly support. It is below the 20%
  breadth sensitivity cutoff. There is no archaeal MAG; the MAGSP0026
  Perkinsela-related reassignment remains valid.
- Associate-virus calls total 631, across 127 libraries and 392 targets.
  Relative to the old cohort after ID translation, 13 old calls disappear
  on six removed targets, 15 calls are gained in original libraries and 99
  occur in Bargmannia. Of the latter, 86 occur in formalin-preserved
  B. elongata T676-SS3. Keep that specimen context when discussing the increase.
- T676-SS3 is the only library whose predefined SSU background groups exceed
  the locked 37-pair read threshold: 288 insect and 71 plant pairs at ≥97%
  identity. Primary grades are retained, with a specimen flag and a descriptive
  cross-domain summary omitting this library in
  `figures/expanded_review/background_sensitivity.tsv`. Explain this limitation
  alongside its many viral calls; the SSU evidence does not by itself invalidate
  the viral matches. The supplied Mańko paper describes a distinct extraction
  protocol for the formalin-preserved nectophore, without identifying the cause
  of the background. Do not claim that the old maximum of 36 background pairs
  holds across the expanded cohort, or retune the threshold after this result.
- The 213-library fish screen totals 9,313,145,671 filtered pairs. Strong
  results remain 32 distinct paired sequences in eight Physalia libraries.
  The eight Bargmannia screens add no strong fish matches. The selected inputs
  are 163 host-depleted libraries (including both Halistemma) and 50 full-read
  libraries. Halistemma reruns use all 131,703/2,352,421 retained pairs and yield
  no qualifying fish matches; they replace the earlier full-read entries.
- Flowcell-based inference includes 144 libraries, of which 122 are Physalia.
  The primary Physalia community test uses the 91 nonempty communities:
  geographic-region R² = 0.1890813593, P = 0.0001; three of five tested
  individual targets have q < 0.05. No target is flagged by the current
  flowcell test. This does not establish absence of contamination.
- Eukaryotic results total 64 nonredundant detections in 40 libraries and
  17 groups, including 42 high-confidence detections. All 61 original calls
  and their grades remain. Three read-only additions occur in T676-SS3:
  fungi, Bacillariophytina and the perkinsid group A31. Omitting the flagged
  library gives 61 eukaryotic, 567 microbial and 545 viral detections. Its
  86 viral targets are not detected elsewhere at the accepted threshold.
  Include the specimen qualification in the eukaryote figure legend and viral
  results rather than presenting these additions as unqualified evidence.
- The refreshed screen comparison has 505 sylph-nominated genomes and 977
  nomination/library pairs. Of 926 pairs represented in the final catalog,
  283 pass presence criteria and 643 do not. Of 568 accepted microbial calls,
  283 were nominated in that library, 29 only in other libraries, and 256 are
  from 24 MAG-only targets. All 77 libraries with the seven focal genomes have
  broad Mycoplasma-related Bracken and phyloFlash support. Use generated
  `figures/screening/summary.json` values and their precise counting units.

## Manuscript files requiring reconciliation

The manuscript's `scripts/sync_assets.py` still hardcodes
`data/results/phase6_analysis/v2` and the old `figures/catalog_review/mag_catalog.tsv`.
Update the selection to the accepted paths in
[reporting_sources.json](../config/reporting_sources.json) before transferring
assets. Retain the historical domain-audit evidence separately from current
catalog counts. The old `domainRetainedMags`/cluster macros must not silently
supply the v3 totals. Transfer provenance should identify these selected sources.

| Manuscript location | Required review |
| --- | --- |
| `sections/abstract.tex`, `introduction.tex`, `data-summary.tex` | Four-study/213-library scope and Mańko et al. citation; current metadata and catalog links |
| `sections/methods.tex`, `supplement/methods.tex` | Halistemma reference/depletion, universal assembly attempts, v3 curation and mapping, expanded host SSU competition, dated cached verification, current metadata coverage |
| `sections/results.tex`, `supplement/results.tex` | Current MAG/target/detection totals, Tyloplasma, limited archaeal reference evidence, changed viral findings, current fish denominator and completed eukaryote/screen comparisons |
| `sections/discussion.tex` | Preserve the false-archaeal lesson while distinguishing the new reference-only signal; incorporate new findings without extending the seven-genome functional scope |
| `supplement/figures-and-tables.tex` and main table/figure legends | All 213 descriptive denominators, new study/hosts, current evidence counts and depth/tissue metadata |
| Figure 1 source | Add the fourth source study and Halistemma route where the diagram specifies cohorts or host references |

The original host-depletion/full-read and capped/uncapped method checks remain
v2 analyses of fixed subsets. Keep their original catalog and denominators
explicit; they were not repeated in this refresh. The user previously ruled out
a costly clean full-sequence rerun. Current checks establish reuse of verified
outputs and successful bounded downstream execution.

Use the completed source tables and baseline comparison to update literal
values as well as generated macros. `figures/overview/study_overview.pdf` is a
retained concept sketch with photo placeholders, not the manuscript's Figure 1;
update the actual manuscript drawing. Rebuild and inspect the
manuscript only in the subsequent authorized manuscript task.
