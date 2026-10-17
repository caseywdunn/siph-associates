# Fish 12S screen: initial assessment

Work on branch `fish-12s-screen`, started 2026-10-17. The intended cohort is all
205 libraries, with Physalia diet as the main biological question. This first
stage tests the proposed regex and compares candidate-read filters. It does not
yet report fish identifications or run the full cohort. Existing sequence
analyses, metadata identifiers, and manuscript results are unchanged.

## Does the original regex cover fishes broadly?

The two alternatives are exact reverse complements. The motif spans 29 bases,
with 23 fixed positions and six wildcard positions:

```text
ACAAAC..GGATT..ATACCCC..TATGC|GCATA..GGGGTAT..AATCC..GTTTGT
```

We tested the September 2026 [MitoFish](https://mitofish.aori.u-tokyo.ac.jp/download/)
snapshot. To avoid treating short amplicons that omit the motif as failures,
the assessment uses annotated 12S sequences at least 700 nt long, extracted from
assemblies at least 14,000 nt long, with only A/C/G/T bases. Identical 12S
sequences within a species are counted once. This leaves 11,516 sequence/species
combinations across 7,236 reference species. Six extractions, one from each major
fish group, were checked against their original GenBank sequences; all matched.

| Fish group | Reference sequences | Original regex | Allow one substitution at a fixed position |
| --- | ---: | ---: | ---: |
| Ray-finned fishes | 10,837 | 10,023 (92.5%) | 10,757 (99.3%) |
| Sharks, rays, and chimaeras | 624 | 0 (0%) | 169 (27.1%) |
| Lampreys | 35 | 34 (97.1%) | 34 (97.1%) |
| Hagfishes | 10 | 0 (0%) | 0 (0%) |
| Coelacanths | 3 | 3 (100%) | 3 (100%) |
| Lungfishes | 7 | 3 (42.9%) | 7 (100%) |

These are reference-sequence coverage estimates, not read-level sensitivity or
an exhaustive survey of fish species. The one-substitution comparison allows no
insertions or deletions. [Detailed coverage](../dev_docs/fish-12s/coverage_summary.tsv)
also records species coverage and finer taxonomic breakdowns.

The original motif is therefore useful for initial candidates but is not a
universal fish screen. A common cartilaginous-fish variant differs at two fixed
positions: `CCAAACTAGGATTAGATACCCTACTATGC`. Even when the motif is present in a
species, reads must span that particular site to match; most reads elsewhere in
12S cannot be recovered with this filter.

Specificity also requires a subsequent identification step. Of five non-fish
mitochondrial controls, *Xenopus laevis* matched the original motif exactly;
human, mouse, cow, and chicken did not. The whole-12S filter described below also
retains the frog control. Neither filter alone identifies a read as fish.

## Performance and read pilot

The read inventory contains 313 paired lane entries representing 205 libraries,
approximately 6.05 TB compressed and 40.1 billion read pairs. We first used the
first million pairs from one lane in each of six libraries, covering all three
studies. WAM-Z97917 and YPM-IZ-111158 were included because earlier SSU analyses
found fish evidence. A fixed lane prefix is a technical pilot, not a
representative subsample for conclusions about absence.

| Method | Median seconds per million pairs | Output |
| --- | ---: | --- |
| Original `zgrep -E`, both files | 3.44 | Matching lines |
| Same command with `LC_ALL=C` | 3.52 | Matching lines |
| SeqKit sequence-only regex, one thread | 16.72 | Matching FASTQ records |
| BBDuk whole-12S reference, four processing threads | 5.06 | Both mates when either matches |

These timings use the same bounded, recompressed inputs and include tool startup.
The methods have different CPU usage; this is not a full-dataset throughput
estimate. Setting `LC_ALL=C` did not improve these runs. SeqKit was slower than
grep here, although it preserves FASTQ records and avoids matching quality lines.

The original regex found no reads in these six million pairs. BBDuk found two
pairs in NA31; both were dominated by simple repeats, not useful fish evidence.
We therefore made a search-only derivative of the 12S reference using
[NCBI DUST](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/dustmasker/README)
(window 64, level 20, linker 1), replacing masked bases with N. It masks 22,575 of
10,960,421 bases (0.206%) while preserving record identities and lengths.
The two repetitive pairs no longer match this derivative.

The candidate filter uses [BBDuk](https://bbmap.org/tools/bbduk), exact 31-mers,
both orientations, no middle-base wildcard, and both mates retained if either
matches. Synthetic controls verify forward and reverse-complement detection,
mate recovery, a cartilaginous-fish variant missed by the regex, and fish reads
outside the motif. Negative reads and quality-only matches are rejected; the
masked reference also rejects poly-G and ATTT repeats. These controls establish
implementation behavior, not sensitivity to fish absent from the reference.

## Recommended next stage

Use the masked whole-12S filter across all libraries, retaining complete read
pairs, names, qualities, source lane, and original analysis sample ID. Keep the
original regex as a comparison applied to the candidates. Before biological
reporting, compare candidate alignments against fish and non-fish references,
resolve ambiguous assignments, inspect duplicates and support outside conserved
sites, and evaluate sequencing-batch contamination. The SSU analysis's 37-pair
threshold has not been calibrated for this marker and should not be transferred.

The main speed improvements are a single read pass, recovery of mates during that
pass, and modest parallelism across lane files. The streaming filters need far
less memory than the pilot shell script's 128 GB. Benchmark a full lane before
choosing production CPU counts and disk concurrency. No full-cohort scan has
been launched during this assessment.

## Reproduce the assessment

Download URLs, sizes, and SHA-256 hashes are in
[reference_downloads.tsv](../dev_docs/fish-12s/reference_downloads.tsv).
Place those files in `data/results/fish_12s/reference/2026-10-17/` and extract
`tables.tar` there. Reference results and large FASTA/FASTQ files are excluded
from Git. Small evidence tables are in [dev_docs/fish-12s](../dev_docs/fish-12s/).

The assessment needs Python, NumPy, and PyArrow. On McCleary:

```bash
module load Arrow/16.1.0-gfbf-2022b
python scripts/fish_12s/assess_reference.py \
  --reference-dir data/results/fish_12s/reference/2026-10-17 \
  --output-dir data/results/fish_12s/assessment/NEW_RUN
python -m pytest -q tests/test_fish_12s_reference.py
```

The masking, controls, and read pilot use Biopython and the existing
`/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python` environment.
The paired-read pilot used BBTools 39.01 from `phyloflash_env` and SeqKit v2.8.0
(the cluster module is named `SeqKit/2.8.1`). The SLURM launcher establishes both
tool paths, including Java. Its output directory is unique per sample and fails
if already present; choose a new output-root for reruns.

```bash
mkdir -p logs/fish_12s
sbatch scripts/fish_12s/pilot.sbatch
sbatch --array=0-1%2 scripts/fish_12s/pilot.sbatch config/fish_12s_positive_pilot.tsv
```

These are the actual initial jobs: `13498984` and `13499007`, all completed.
They use one million pairs per library. Full commands, source/input hashes, tool
versions, timings, and output hashes are retained under
`data/results/fish_12s/{assessment,pilot,controls,masked_reference}/` in `run.json`
and `commands.json` files. Original script snapshots are in
`data/results/fish_12s/source_snapshots/initial/`. The direct masked recheck of the
two candidate pairs is recorded in `masked_pilot/2026-10-17/`; an initial missing
Java PATH failure is retained separately from the successful retry.

The reference assessment and standalone pilots have no dependencies on the
existing Snakemake targets and cannot schedule earlier sequence analyses.
