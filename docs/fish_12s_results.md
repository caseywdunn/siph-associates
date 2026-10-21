# Fish 12S screen across the expanded cohort

All **213 libraries** completed the screen: all retained host-depleted reads for
151 Physalia and ten Nanomia septata libraries, and up to 200 million raw pairs
before quality filtering for each of the other 52 libraries. The latter include
Erenna, Rhizophysa, both Halistemma and the eight new Bargmannia libraries.
Halistemma fish screening preceded host depletion and used the full read route.

The [generated summary](../figures/fish_12s/summary.json) and
[library table](../figures/fish_12s/sample_summary.tsv) combine the original
cohort and completed Bargmannia reports without repeating searches. The original
205 libraries contributed 8,146,578,446 screened pairs; Bargmannia adds
1,430,751,915, for **9,577,330,361 pairs**. Full-read libraries contributed
9,160,588,614 quality-filtered pairs.

The 32 long, close competitive fish matches remain confined to eight Physalia
libraries. The eight Bargmannia libraries add no qualifying pairs, consistent
with the absence of qualifying pairs in the other original full-read libraries
and ten N. septata libraries. There are 81 weaker or ambiguous fish-similar pairs
across the cohort (64 original and 17 Bargmannia); these are not strong matches.

The [positive-pair table](../figures/fish_12s/positive_pairs.tsv) retains both
mates' alignment evidence, near-best reference names, recruitment methods,
flowcells and exact paired-sequence hashes. All 32 paired sequences are distinct.
This checks exact duplication; it does not establish independent original DNA
molecules, exclude sequencing contamination or count prey individuals. Closest
reference names remain identification hypotheses, especially with near-tied hits.

BBDuk recovered all 32 pairs, Bowtie2 and SortMeRNA each recovered 31, and the
original regex recovered seven. Strong matches require at least 100 nt at 95%
identity, fish-only near-best hits within 5% of the best bit score, an adequately
resolved hit-list tail, and no strong non-fish competing mate. RefSeq membership
uses accession-linked GenBank annotations. This reproduces the earlier result
without relying on fish genus names. Whole mitochondrial references can extend
alignments beyond 12S. The SSU analysis's 37-pair threshold is not applied here.

## Reproduce reporting

The [execution guide](prey_marker_screen.md) describes the searches and controls.
From the repository root, with the existing Snakemake Python environment:

```bash
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python \
  scripts/fish_12s/finalize_report.py
```

This reclassifies preserved host-depleted alignments, verifies report hashes,
checks cohort membership and pair counts, and combines the 44 full-read reports.
It performs no new sequence searches. Complete candidate evidence is retained at
`data/results/fish_12s/final_review/2026-10-19/all_candidate_pairs.tsv`;
the small reporting bundle and provenance are in `figures/fish_12s/`.
Failure to recover a qualifying read does not establish that a host never eats fish.
