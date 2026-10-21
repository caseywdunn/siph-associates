# Fish 12S screen across the cohort

All **213 libraries** completed the screen: all retained host-depleted reads for
151 Physalia, ten Nanomia septata and two Halistemma rubrum libraries, and up to
200 million raw pairs before quality filtering for each of the other 50 libraries.
All 163 host-depleted libraries use the same search and classification methods.

The [generated summary](../figures/fish_12s/summary.json) and
[library table](../figures/fish_12s/sample_summary.tsv) combine the selected
per-library reports: **9,313,145,671 pairs** in total, comprising 419,225,871
host-depleted pairs and 8,893,919,800 quality-filtered full-read pairs.

The 32 long, close competitive fish matches remain confined to eight Physalia
libraries. No other library has a qualifying pair. There are 81 weaker or
ambiguous fish-similar pairs across the cohort; these are not strong matches.

Halistemma CWD12 and CWD17 contributed 131,703 and 2,352,421 retained pairs,
respectively. CWD12 yielded no candidates; CWD17 yielded three candidates, none
classified as fish similarity. Their prior full-read screens are excluded from
the combined report. The [execution record](../dev_docs/fish-12s/2026-10-21-halistemma-host-depleted.md)
preserves the comparison and validation.

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

This verifies report hashes, checks cohort membership and pair counts, and
combines 163 host-depleted and 50 full-read library reports. It selects the two
Halistemma depleted reports in place of their full-read reports and performs no
new sequence searches. Complete selected candidate evidence is retained at
`data/results/post_expansion/2026-10-21/reporting/fish_12s/all_candidate_pairs.tsv`;
the small reporting bundle and provenance are in `figures/fish_12s/`.
Failure to recover a qualifying read does not establish that a host never eats fish.
