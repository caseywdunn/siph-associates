# C. acnes and human-compatible prefix review

## Scope and design

Coauthor review requested a human contamination screen. Once the 213-library
bounded screen was complete, compare it descriptively to accepted BAC00019
*C. acnes* detections. Preserve selected-prefix denominators, reporting grades
and uncertainty about both source and taxonomic specificity. No life-stage
analysis, new read processing or remapping is performed by these entry points.

The four largest human-compatible counts occur in specimens with accepted
*C. acnes* detections. Median counts are 61 versus 3 pairs per million selected
raw pairs in specimens with versus without an accepted *C. acnes* detection.
Omitting the independently flagged formalin specimen T676-SS3 leaves a median
of 50 for the *C. acnes* group. NA37 has accepted *C. acnes* but no retained
human-compatible pairs. No formal test, contamination cutoff or taxonomic
abundance correction is applied. The nt representative-sequence review remains
separate.

A single targeted community sensitivity was then requested: omit BAC00019 from
the two estimable primary community contrasts, keeping sample definitions and
all statistical settings. First reconstruct both baseline results exactly;
then remove only that target. This does not overwrite original results or
remove specimens based on human screening.

## Execution

Both workflows ran within the existing one-core compute allocation 13489566
on 2026-10-22. The descriptive join completed in under one second, and the
community baseline reconstruction plus omission sensitivity took twenty seconds.

- Join log: `.snakemake/log/2026-10-22T143406.863266.snakemake.log`.
- Sensitivity log: `.snakemake/log/2026-10-22T143630.060434.snakemake.log`.
- Reports: `figures/human_cutibacterium/` and
  `figures/cutibacterium_sensitivity/`.
- Reproduction and interpretation:
  [descriptive comparison](../docs/human_cutibacterium.md) and
  [community sensitivity](../docs/cutibacterium_sensitivity.md).

## Validation and result

Four focused Python tests protect one-to-one specimen joins, accepted-only
bacterial evidence, human raw/QC denominators, impossible counts and flagged
specimen omission. The join reconciles exactly to the human-screen cohort:
213 million selected raw pairs, 194,926,432 QC pairs and 4,563 compatible pairs
across 178 libraries. Sixteen specimens have accepted *C. acnes* detections.

Both primary community baselines reproduce their original counts, model
degrees of freedom, R-squared, pseudo-F and P values within 1e-10. After omission,
membership is unchanged. The southwestern-Pacific host comparison retains
R-squared 13.3%, P 0.0001 and q 0.0002. The Pacific *P. utriculus* regional
comparison has R-squared 7.1% and P = q = 0.1539. No R warnings occurred.

Python formatting/lint, provenance hashes and output-source snapshots were
checked. Both rule graphs are registered for the existing graph checks. The
new R script uses the retained statistics environment; lintr and styler are
unavailable there. The exact baseline reconstruction is the main numerical
regression check.
