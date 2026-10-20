# Fish 12S screen across the original cohort

All 205 libraries completed the screen: all retained host-depleted reads for
151 Physalia and ten Nanomia septata libraries, and up to 200 million raw pairs
before quality filtering for each of the other 44 libraries. The latter include
Erenna and Rhizophysa. Array 13512346 completed all 44 tasks successfully.

The [generated summary](../figures/fish_12s/summary.json) and
[library table](../figures/fish_12s/sample_summary.tsv) report 8,146,578,446
screened pairs. There are 32 long, close competitive fish matches in eight
Physalia libraries, with no additional qualifying pairs in the 44 full-read
libraries or ten N. septata libraries. The full-read libraries contributed
7,729,836,699 quality-filtered pairs. Sixty-four other pairs across the cohort
have weaker or ambiguous fish similarity and are not strong matches.

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
