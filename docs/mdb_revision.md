# Summaries for Martina Dal Bello's manuscript comments

This reporting step reuses accepted v2 detections and the completed gene-content
comparison. Run from the analysis repository with Python 3.11 or later
(standard library only):

```sh
python3 scripts/summarize_mdb_revision.py
```

Outputs in `figures/mdb_revision/` include specimen-level joins, descriptive
co-occurrence tables, selected nutritional/defense annotations, generated LaTeX
tables and phage text, and an input/output SHA-256 provenance record. The
manuscript's explicit asset sync copies these reviewed outputs. No workflow
rule, sequence result, detection threshold or genome catalog changes.

## Questions and scope

- **Co-occurrence:** one specimen per library; retain all 205, including zero
  detections. The primary descriptive comparison is within Physalia. Count a
  specimen once per exact GTDB genus (Vibrio, Pseudoalteromonas, Alteromonas,
  Photobacterium), and once for their union. These are candidate ecological
  groups, not experimentally or genomically verified chitin-degrading strains.
  Reuse Metamycoplasmataceae and clade indicators from the collection-context
  table; independently verify the family and Vibrio indicators against the
  detection table. Separate clade A, clade B and DT-68 summaries are retained.
- **Confounding:** provide single-flowcell, tentacle and P. utriculus tentacle
  subsets, plus exact host-species/region/tissue/flowcell strata in eligible
  Physalia. Retain input-depth medians. These are descriptive diagnostics,
  not covariate-adjusted estimates or significance tests. Genus detection
  does not measure abundance enrichment or demonstrate nutritional exchange.
- **Nutritional limitations:** extract stored amino-acid biosynthetic modules
  under both annotation tiers, excluding plant melatonin and ethylene production
  (which consumes methionine). Do not reinterpret a missing module row as a
  measured zero. Display four precursor-dependent routes with returned focal
  assignments (cysteine, methionine, histidine, arginine); keep the full selected
  module set in TSV. Do not add overlapping modules to count auxotrophies or
  infer which organism supplies required nutrients.
- **Defense:** display the prespecified Cas and restriction–modification KOs,
  preserving tier differences. Count assigned genes, not complete systems.
  A comprehensive defense survey remains a separate potential annotation task.
- **Genome quality:** the five near-complete genomes support cautious absence
  discussion. MAGSP0029 and MAGSP0031 remain presence-only in supporting TSVs.
- **Phages:** generate the existing six-pair result directly from the retained
  host-link table, including the distinction between spacer linkage and
  contemporaneous co-detection.

These are retrospective questions prompted by co-author review. No test was
selected for significance. Physalia's fish diet precludes presenting direct
crustacean feeding as its general source of chitin; the manuscript treats
substrate availability and bacterial localization as unresolved.
