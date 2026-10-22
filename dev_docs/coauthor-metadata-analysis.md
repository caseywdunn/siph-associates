# Coauthor metadata analysis: execution and validation

**Scope:** post hoc tissue/depth summaries, matched bacterial community
comparisons, fish-marker tissue context and host-reference summaries.
**State:** validation passed, 2026-10-22. No sequence processing, catalog
reconstruction, mapping or life-stage analysis was run.

The design and reproduction commands are in
[docs/coauthor_metadata_analysis.md](../docs/coauthor_metadata_analysis.md).
The final selected execution is
[`13489566-20261022T172724Z`](../figures/coauthor_metadata/provenance.json),
run with one core in the existing agent compute allocation. Snakemake completed
all four steps in about 33 seconds; its execution log is
`.snakemake/log/2026-10-22T132651.497755.snakemake.log`. The first trial retained
the preceding viral classification and NA37 identification and is superseded
by this final refreshed execution. No statistical choices were selected by
their significance. The final reporting repair filled the residual-df output
column, and the dispersion diagnostic uses square-root Jaccard to avoid
negative squared distances; the raw-Jaccard PERMANOVA was unchanged.

## Result and interpretation

Within southwestern-Pacific tentacle samples, host species accounted for
13.08% of detected bacterial-profile variation after adjustment for mapping
input pairs and physical flowcell (30 nonempty profiles, permutation
P=0.0001, BH q=0.0002). The association persisted at 5% and 20% breadth
(q=0.0004 and 0.0466). The dispersion check did not identify differing spread
(P=0.9995). This is an observational host association, with collection-date
and finer-scale geographic differences still possible.

Within *P. utriculus* tentacles, the three Pacific regions accounted for
6.77% of variation in 32 nonempty profiles, but the primary test was not
significant (P=q=0.2016). The result varied with the breadth criterion:
P=0.0092 at 5% and P=0.1371 at 20%. The broader western-Indian comparison
was rank-deficient after empty-profile exclusions. Therefore the original
all-Physalia geographic result should not be described as evidence that
geography is independently established within this host species.

Fish 12S matches occurred in two of three zooid samples and six of 142
tentacle samples. The zooids accounted for 19 of the 32 distinct matching
read pairs. These small, uneven samples support a descriptive tissue
observation only. Other tissue categories had no retained fish 12S matches.

All 213 depths are recorded or explicitly derived, but the major DT-68
contrast in *N. septata* remains entirely coupled with locality: seven deeper
California specimens versus three shallow Washington specimens. The most
closely matched tissue comparison has only twelve distinct blocked label
allocations before empty-profile exclusion. No independent depth or tissue
effect was claimed.

Current host-depletion medians were 98.53% in Physalia, 99.90% in *N. septata*
and 99.23% in *H. rubrum*. The historical ten-library fixed-v2 comparison
retained 21 bacterial detection pairs in both read treatments, one only
after depletion and three only with full reads. Its combined viral catalog
retained 22 pairs in both and 31 only with full reads; those counts include
historical endogenous candidates and are not current associate-virus totals.

## Checks completed

- Five focused tests passed for one-to-one mapping joins, duplicate rejection,
  marker/genome denominator separation, blocked tissue allocations and
  depth/locality confounding.
- Ruff format and lint passed for all new Python sources/tests; the R source
  parsed and completed without warnings. R `lintr`/`styler` are not installed
  in the retained statistics environment; no claim of those checks is made.
- The full cohort contained 213 unique libraries, NA37 was *N. bijuga*,
  and the accepted viral total was 624 after the reviewed classification
  correction. The fish total was 32 paired sequences.
- Every input, source and output hash in the final provenance record matched
  the current artifact. Actual R/package versions are in `session_info.txt`;
  source snapshots cover the uncommitted implementation.
- All eligible/included/excluded community counts reconciled; tested contrasts
  retained residual degrees of freedom. Primary and both threshold
  sensitivities completed.
- Current host-depletion statistics matched their original completion/hash
  record and all 163 reference-mapped libraries. Historical paired counts
  reconciled separately in both domains and all ten libraries.
- Both figures were visually inspected; axis labels, legends and fractions
  were readable. Generated values provide the manuscript numbers.
- The bounded workflow dry run and generated rule-graph freshness check
  completed. The entry point has no sequence-processing rules.

The scientific limitations are those in the design document. Results describe
detected targets and cannot convert non-detection into demonstrated biological
absence. The reviewed evidence and read-handling procedures are retained.
