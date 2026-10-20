# Fish 12S screen: initial assessment

Work on branch `fish-12s-screen`, started 2026-10-17. The intended cohort is all
205 libraries, with Physalia diet as the main biological question. The 161-library
host-depleted screen is complete, and the approved expansion covers the other 44
libraries with up to 200 million raw pairs each. This page preserves the initial
method assessments; the [current execution guide](prey_marker_screen.md) documents
the cohort workflows, controls and preliminary result tables. Earlier raw-lane
pilots below are superseded by the complete retained host-depleted inputs.

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

## Positive-control assessment

Before expanding to all libraries, test recovery in libraries with previous fish
evidence. Casey's pilot regex recovered fish 12S from YPM-IZ-110906 and
YPM-IZ-110907. The comparison set includes these two and three other libraries
with retained fish SSU detections: WAM-Z97917, YPM-IZ-111156, and YPM-IZ-111158.
YPM-IZ-110907 also has retained fish SSU evidence (131 read pairs and an assembled
SSU); YPM-IZ-110906 has only a trace vertebrate SSU signal (20 pairs). These are
selection criteria, not new 12S results. Nuclear SSU evidence does not guarantee
that mitochondrial 12S reads will be recovered.

[The manifest](../config/fish_12s_known_positive.tsv) lists every lane for these
five libraries: nine paired entries, 98.03 GB compressed. The first assessment
uses complete lanes, with independent sequence-only searches for the original
regex and exact 25-mers from the masked whole-12S reference. Both mates are kept
for BBDuk candidates. Exact 31-mer membership is then assessed on the 25-mer
candidates: every exact 31-mer necessarily contains exact 25-mers, so this gives
the same retained pairs as scanning the original inputs with k=31. Its timing
measures only the candidate recheck, not full-input k=31 throughput. Regex-only
reads retain their headers and qualities; their unmatched mates are not included
unless BBDuk recovered the pair.

Bowtie2 (`--very-sensitive-local`) and SortMeRNA independently search the complete
L003 lane of each of YPM-IZ-110906 and YPM-IZ-110907 against the same masked 12S
reference. They receive raw inputs, not BBDuk-selected reads. SortMeRNA retains
both mates when either aligns. Bowtie2's SAM records preserve aligned-read
evidence, including alternative-species ambiguity in MAPQ; a single reported
alignment is not a taxonomic assignment. Synthetic checks recover four fish
examples and the homologous frog example with both methods, while rejecting
random, quality-only, and repetitive negative controls.

Candidates from the regex/BBDuk union are compared by BLASTN with unmasked
MitoFish 12S sequences and the complete
[NCBI mitochondrial RefSeq collection](https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/),
file `mitochondrion.1.1.genomic.fna.gz` dated 2026-09-03. This adds mitochondrial
sequences from other organisms. It is not an unrestricted search of all NCBI
nucleotide records, and closest reference matches should not automatically be
reported as prey species. Alignment checks use word size 11, DUST, E-value
1e-8, one HSP per subject, and up to 200 subjects per query. The query-ID table
connects every candidate to its original read, mate, lane, and filter membership.

The full-cohort method remains undecided pending positive-control recovery,
reference-held-out sensitivity tests, and comparisons of computational cost.
Before biological reporting, resolve ambiguous assignments, inspect duplicates
and support outside conserved sites, and evaluate sequencing-batch contamination.
The SSU analysis's 37-pair threshold has not been calibrated for this marker and
should not be transferred.

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

## Run the positive-control assessment

Run from the repository root on McCleary. These commands consume the existing
MitoFish reference and masked derivative described above. They are separate from
the main Snakemake workflow. Select fresh output directories for repeated runs;
the scripts fail if their output directories already exist.

SortMeRNA 7.1.0 was installed in an isolated environment under the ignored
`data/results/fish_12s/software/` directory. Its exact package specification is
[sortmerna_environment.txt](../dev_docs/fish-12s/sortmerna_environment.txt):

```bash
module load miniconda/24.7.1
conda create --prefix data/results/fish_12s/software/sortmerna_env \
  --file dev_docs/fish-12s/sortmerna_environment.txt
```

Download the mitochondrial RefSeq file into
`data/results/fish_12s/reference/refseq_mito_2026-09-03/` using the filename above;
the pilot's download headers and hashes are retained with the run. The directory
listing is mutable, so check the archived input checksum before reproducing this
specific run with a later download. The reference-preparation job uses Bowtie2
2.5.1 and BLAST+ 2.15.0 cluster modules. The read screens use BBTools 39.01 and
SeqKit v2.8.0 as above; launchers set their module and environment paths.

```bash
mkdir -p logs/fish_12s
sbatch scripts/fish_12s/prepare_positive_references.sbatch
sbatch scripts/fish_12s/screen_positive.sbatch
# After reference preparation completes:
sbatch scripts/fish_12s/compare_aligners.sbatch
# After both reference preparation and the read screens complete:
sbatch scripts/fish_12s/align_positive_candidates.sbatch
```

To run the same implementation on a small control, supply a manifest with paired
fixture paths and their byte sizes, set `--array=0`, and pass a fresh output root
as the second script argument. See
[the execution record](../dev_docs/fish-12s/positive_pilot.md) for the actual jobs,
dependencies, checks, and output paths. Candidate FASTQ files, alignments, source
snapshots, and execution logs are retained under `data/results/fish_12s/` and are
excluded from Git. This pilot changes no sample IDs or existing analysis targets.
