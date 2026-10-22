# Associate analysis review, 2026-10-22

## Scope

This record covers analysis changes to the 213-library survey. Manuscript
sources, drafts and detailed coauthor feedback are maintained separately in
`../manuscript_siph_associates`.

Completed sequence assemblies, microbial/viral mapping targets, raw mappings
and stable analysis sample identifiers were preserved. No full-cohort
microbial/viral remapping was launched. All six frozen sequence and mapping
input files retained their original checksums. Recorded Physalia life stages
are retained in source metadata but excluded from analyses because their
assignment was inconsistent.

## Metadata and reporting

- Author review identifies NA37 (`Ahuja2026__NA37`, `YPM:IZ:111590`) as
  *Nanomia bijuga*. NA29 remains unresolved. The publication-name override
  records the correction without changing frozen sequence inputs.
- Three Negarnaviricota targets are now possible endogenous candidates.
  The reviewed viral reporting catalogs contain 404 associate and 160
  endogenous-candidate groups; associate detections total 624 in 122 specimens.
  Original mapping measurements and evidence grades are unchanged.
- `workflow/coauthor_reporting.smk` rebuilds reviewed reporting without
  sequence-analysis producers. Its validation passed with zero errors.

## Tissue, depth and host/region comparisons

`workflow/coauthor_metadata.smk` uses completed detection tables and updated
metadata. It checks confounding and restricts targeted comparisons to tentacle
samples with eligible flowcell records, controlling for read effort.

- The southwestern-Pacific host comparison is supported at all three detection
  thresholds (primary R² 13.1%, P=0.0001).
- The within-*P. utriculus* Pacific regional comparison is threshold-sensitive
  (primary R² 6.8%, P=0.2016). Indian-Ocean contrasts are confounded with flowcell.
- Tissue and depth sampling does not support general independent effects.
  Fish 12S occurs in 2/3 zooid and 6/142 tentacle samples; zooids supply 19/32
  qualifying pairs. These are not abundance estimates.
- In *N. septata*, DT-68 occurs in seven California specimens collected at
  217–662.29 m and none of three Washington specimens from 0–20 m.
  Geography and depth cannot be separated.
- Host-reference mapping removes median 98.53% of pairs in Physalia,
  99.90% in *N. septata* and 99.23% in *H. rubrum*. Historical fixed-catalog
  comparisons do not isolate a causal effect on genome discovery.

## Human-read assessment and contamination sensitivity

The independent human screen selects one million raw pairs per library as
equal per-source-file prefixes. Its 213-library mapping and subsequent nt
review completed in jobs 13734795 and 13735568.

- Of 483,068 concordant hg38 pairs, 4,563 passed paired full-length identity,
  mapping-quality and sequence-complexity filters.
- Controls recovered 89.87% of human QC pairs and none of 29,355 simulated
  siphonophore QC pairs.
- After within-library deduplication and a 500-pair/library cap, 3,384 pairs
  were compared with nt. Both mates of 2,060 pairs had human-only qualifying
  near-best matches; 1,241 retained alternatives, 29 had incomplete evidence
  and 53 were unresolved. One pair matched only a synthetic-construct label.
- Fractions describe the selected prefixes, not complete-library abundance.
  Public outputs contain aggregate summaries rather than candidate sequences.
- The 16 *Cutibacterium acnes*-positive specimens had median 61 human-compatible
  pairs per million selected pairs versus three in other specimens. This
  descriptive association supports a contamination concern without identifying
  its source. Excluding *C. acnes* preserves the southwestern-Pacific host
  comparison (R² 13.3%, P=0.0001); the primary within-*P. utriculus* comparison
  remains unsupported (P=0.1539).

## Validation and reproducibility

Source/input/output checksums, focused tests, lint and bounded workflow dry
runs are recorded with each analysis. The 255-file publication asset bundle
passed its hash checks and was transferred to the separate manuscript
repository. No sequence analyses were rerun for that transfer.

Close-relative deletion tests, simulated detection curves and microbial SNP
analyses were not performed. Source-library recovery and distant decoys do not
calibrate close-relative misassignment or low-abundance sensitivity.

Detailed methods and execution records:

- [Metadata methods](../docs/coauthor_metadata_analysis.md) and
  [execution](coauthor-metadata-analysis.md).
- [Viral reporting methods](../docs/viral_reporting_review.md) and
  [execution](viral-reporting-review-2026-10-22.md).
- [Human-read methods](../docs/human_read_screen.md) and
  [validation](human-screen/2026-10-22-validation.md).
- [Human/C. acnes comparison](../docs/human_cutibacterium.md),
  [community sensitivity](../docs/cutibacterium_sensitivity.md) and
  [execution](cutibacterium-human-review.md).
