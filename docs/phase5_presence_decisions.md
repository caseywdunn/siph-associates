# Phase-5 presence rules

Approved by Casey Dunn on 2026-09-29, after cohort mapping against catalog v1
and before any presence was interpreted. The rules are frozen in
`config/phase5_presence.json`. The control tables are regenerated from
accepted mapping outputs with `python3 scripts/phase5_presence_controls.py`,
which writes to `data/results/phase5_mapping/cohort/presence_controls/`.

The mapping filters were fixed before mapping (`config/phase5_mapping.json`):
bwa-mem `-M -T60`; primary, properly paired alignments at MAPQ ≥30;
duplicates removed; CoverM counts reads at ≥95% identity over ≥75% of the
read. Breadth is the fraction of a genome or vOTU covered by such reads.

## 1. Bacterial and archaeal presence

**Decision.**
- **Validated:** ≥10% genome breadth and ≥100 reads.
- **Trace or ambiguous:** 1–10% breadth. Reported, never claimed as presence.
- Phase 6 repeats the headline analyses at 5% and 20% breadth.

**Controls** (`rule_grid.tsv`, `positives.tsv`, `decoy_recruitment.tsv`):
- **Positives.** Each of the 38 Phase-3 MAG species, in the library it was
  assembled from. Where a MAG species was represented in the catalog by a
  reference genome, the positive is that genome.
- **Negatives.** The 20 decoy genomes in all 205 libraries (4,100 pairs).

| Breadth ≥ | Positives passing | Decoy–library pairs passing |
|---|---|---|
| 1% | 38/38 | 1/4,100 |
| 5%, 10%, 20% | 38/38 | 0/4,100 |
| 50% | 37/38 | 0/4,100 |

Positive breadth ranges from 41.9% to 100%. The highest decoy breadth is 2.2%
(DEC018 in Ahuja2026 NA31, 436 reads). All other decoy–library pairs stay at
or below 0.25%, even with up to 37,411 reads: decoy reads pile onto a few
conserved loci. Read minimums of 10–100 do not change any result, and every
call at ≥10% breadth has more than 100 reads.

**Reasoning.**
- The plan's starting rule (D0.2) is supported by complete separation, with
  wide margins on both sides. Positives sit fourfold above the threshold and
  the worst decoy fivefold below it.
- The 100-read minimum costs nothing here and protects future catalog versions
  against a short genome reaching 10% breadth from a handful of reads.
- These controls establish specificity against unrelated lineages. Reagent and
  skin contaminants are real genomes that genuinely map (for example
  *Cutibacterium acnes* at 84% breadth in Ahuja2026 NA29). They are graded
  *probable contaminant* from breadth-by-batch patterns in Phase 6, not by the
  presence rule.

**Result** (`bacterial_breadth_bands.tsv`): 795 validated genome–library pairs
(193 genomes, 161 libraries) and 1,433 trace pairs, out of 93,685
(457 genomes × 205 libraries).

## 2. Divergent-strain annotation

**Decision.** Report the breadth ratio for every call: observed breadth divided
by the breadth expected for uniform coverage at the read-derived depth,
1 − e^(−depth), where depth is reads × 150 / genome length. Calls with a ratio
below 0.5 are annotated *divergent-strain match*. The ratio is not a presence
filter.

**Evidence** (`breadth_ratio_calls.tsv`): the ratio is continuous among the 795
calls, with no natural break. 254 calls have a ratio below 0.5, and 236 of
those are MAG species, mostly novel Mycoplasmatales in *Physalia* libraries
other than their source library. For example, BAC00031 reaches about 11%
breadth at 2–9× read depth in several libraries. At that depth a genome
actually present would be nearly fully covered.

**Reasoning.**
- A low ratio means reads concentrate on part of the genome. This is what
  happens when a related strain shares only some regions at ≥95% identity, or
  when conserved loci recruit reads.
- The lineage is detected, but the catalog genome is not the strain present.
  That is scientifically meaningful and should be reported, not discarded.
- As a filter it would remove a true positive (MAGSP0031, ratio 0.42, in its
  own source library) without a principled cutoff.
- Read-derived depth is used because CoverM reports a mean depth of 0 for
  targets whose reads pile onto one locus, which would hide exactly this
  pattern.

## 3. Viral presence

**Decision.**
- **Present:** ≥75% vOTU breadth.
- **Partial or related-virus match:** 10–75% breadth. Reported, not claimed as
  presence.

**Evidence** (`viral_source_summary.tsv`, `viral_source_breadth.tsv`):
- **Positives.** Associate vOTUs in the library their representative was
  assembled from. There are no viral decoys.
- Median breadth is 100%, and 300 of 319 vOTUs (94%) reach ≥75%.
- 18 of the 19 below 75% share 33–98% of their length, at ≥98% identity, with
  another vOTU representative. Reads in the shared regions map ambiguously
  (MAPQ <30) and do not count toward breadth.
- Across the cohort, 530 vOTU–library pairs reach ≥75% breadth and 880 fall
  in the 10–75% band.

**Reasoning.**
- ≥75% breadth is a widely used vOTU detection threshold, and it recovers 94%
  of source-library positives.
- The failures are a known cost of competitive mapping among phages that share
  near-identical modules, not a threshold problem. Lowering the threshold
  would admit partial matches from related viruses, which this catalog
  demonstrably contains.
- This is documented as a sensitivity limit: related phages sharing modules
  can be under-called.
- One such pair (VOTU00053 / VOTU00025) is borderline under the MIUViG
  clustering rule. CheckV's `aniclust` evaluates it from the longer sequence,
  where identity is 94.87%, just below 95%, so the two remain separate vOTUs.
