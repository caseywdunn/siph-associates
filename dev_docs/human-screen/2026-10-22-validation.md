# Human-read audit validation, 2026-10-22

The analysis is independent of the completed microbial and viral mappings. It
uses original, pre-depletion reads to avoid confusing enrichment in retained
reads with the proportion of human-compatible reads in a library.

## Completed control run

SLURM job `13732419` completed the four known-origin controls. The current
settings recovered 89.87% of 9,991 human pairs surviving QC. None of 29,355 QC
pairs simulated from Physalia, Nanomia septata or Halistemma rubrum survived the
human alignment filters. The simulated reads contain one substitution per mate.
The individual counts are retained in [control_summary.tsv](control_summary.tsv).
The initial gate requires at least 70% human recovery and no siphonophore false
positives; it is a method check, not a universal sensitivity estimate.

The complete run outputs are in
`data/results/human_screen/2026-10-22-v2/`. Executed source, versions and exit
status are in `data/results/human_screen/executions/13732419/`.

## Rejected pilot settings

Job `13731699` stopped before processing reads because loading fastp after
Bowtie2 replaced the GCC runtime required by Bowtie2. The launcher loads the
older fastp module first and then Bowtie2; both tools run successfully.

Job `13731866` tested a deliberately restrictive Bowtie2 score minimum
(`L,0,-0.12`) with MAPQ 30 and the explicit identity/length filters. Although all
9,991 human control pairs aligned, only one passed the final filter. Bowtie2's
MAPQ depends on its allowed score range, so tightening the range made the MAPQ
filter unexpectedly reject near-exact, uniquely placed reads. Those biological
sample counts must not be used.

The validated configuration retains Bowtie2's normal end-to-end score range
(`L,0,-0.6`) and applies the independently specified ≥99% identity, full-length,
proper-pair and MAPQ 30 filters afterward. The rejected run remains under
`data/results/human_screen/2026-10-22/`, with its source snapshot and logs.

Twelve unit checks cover the pair decision, including divergence, clipping,
short matches, missing mates, low MAPQ, low complexity and improper pairing.
Four of the checks specifically prevent truncated searches, unknown taxonomy,
near-best nonhuman alternatives, or conflicting mates from being silently
classified as uniquely human. Ruff formatting and lint checks pass. The full workflow is separately scoped;
its graph cannot reach any associate catalog or mapping rule.

## Interpretation boundaries

One million pairs per library are selected as equal per-file prefixes, not a
uniform random sample. The resulting fractions apply only to these prefixes.
The screen does not estimate a complete library's human abundance, identify
individual humans, or automatically remove any reads. Exact sequence-pair
counts and one-megabase reference windows help identify duplication and highly
localized conserved matches.

Candidate confirmation searches representative distinct pairs against the local
2023 nt database and retains near-best nonhuman alternatives. Results remain
human-compatible until that review. An upper bound of 500 distinct candidate
pairs per library applies only to the taxonomic validation, not the mapping
count. No taxonomic claim should treat missing matches in the 2023 database as
proof that newer references would fail to match.

## Comparison with the supplied workflow

The authenticated GitHub API retrieved the exact user-supplied
`dunnlab/illumina_qc/main/Snakefile`. Its SHA-256 is
`945392e92b3c7fbafdbcab04ed3158b5464311d4843aaad92344f2eb8065998f`;
the source is preserved at
`data/results/human_screen/method_sources/illumina_qc.Snakefile`.
`map_human` runs BWA MEM, then requires MAPQ 30, primary proper-pair alignments
and absolute template length 200–1,000. It counts passing mate alignments.
The rule links Casey's
[mapper-stringency note](https://gist.github.com/caseywdunn/e401f94809984bb7e505450b24b7a9be).

This implementation uses Bowtie2 end-to-end mapping with a fixed maximum insert
size, disallows discordant and mixed matches, and requires the complete filter
on both mates of each counted pair. It does not impose the 200-base minimum
insert because genuine fragmented DNA, especially in degraded material, can
produce shorter proper pairs. Full-length, near-exact identity and specificity
review are the safeguards for those pairs. Per-file prefixes represent every
source file but are not a uniform sample of all reads.

## Completed biological pilot

Job `13733275` completed the eight-library pilot with the validated settings.
T676-SS3 has 937 human-compatible pairs among one million selected raw pairs
(0.0937%), representing 681 distinct sequence pairs across 23 reference sequences
and 576 one-megabase windows. The other seven pilot libraries have 13 such pairs
combined (0–6 each). These are mapping results pending broad nt review, not
final human-specific contamination estimates. Full reference/index hashes and
the exact executed source are in the job's execution directory.


## Completed cohort mapping audit

Job `13734795` completed all 213 libraries with 64 allocated cores and 128 GB.
SLURM reports 15 minutes 16 seconds elapsed and 56,381,916 KB maximum RSS.
There were 213,000,000 selected raw pairs and 194,926,432 pairs after QC.
The stringent mapper retained 4,563 human-compatible pairs in 178 libraries.
The largest selected-prefix fractions were NA31 (1,058 pairs; 0.1058%),
T676-SS3 (937; 0.0937%) and D666-D10 (526; 0.0526%). These are proportions
of the selected prefixes only. The audit did not remap any associate catalog.

Pilot nt confirmation `13734768` completed 513 distinct candidate pairs:
322 had both mates matching only human among near-best, nontruncated hits;
185 retained human and nonhuman alternatives, two lacked complete qualifying
hits and four were unresolved because of truncation or unknown taxonomy.
There were no pairs with discordant human-versus-nonhuman best-mate outcomes.
For T676-SS3, the corresponding counts among 500 assessed pairs were 313,
182, one and four. Alternative reference labels include eukaryotic synthetic constructs (taxid
111789), close primates and conserved sequences from other animals. Construct
labels are not evidence for a separate biological organism; they remain
alternatives in this deliberately conservative review.

The shared nt database provides taxonomic IDs but lacks its binary taxdb name
lookup. Readable names are therefore resolved from the local names.dmp after
searching. This reporting fix does not change matches or classifications.

## Completed cohort taxonomic review and report

Job `13735568` completed the full candidate review with 16 allocated cores and
64 GB, running from 14:32:25 to 14:52:45. Of 3,384 assessed candidate pairs,
2,060 had both mates matching only human among qualifying near-best hits,
1,241 retained alternative reference labels, 29 had incomplete qualifying
matches and 53 remained unresolved because of hit truncation or unknown
taxonomy. One pair had only eukaryotic synthetic construct references among
its near-best matches; this label does not establish a nonhuman biological
organism. No pair had discordant human-versus-nonhuman mate classifications.
The review included every distinct candidate pair up to the deterministic
500-pair cap per library. It does not extrapolate human-specific counts to the
remaining pairs or to the full libraries.

SLURM reports 20 minutes 37 seconds elapsed and 67,109,340 KB maximum RSS for
the nt job. BLAST used nearly the entire memory allocation; larger query sets
should request additional memory even though this bounded audit completed.

The initial 483,068 concordant hg38 mappings fell to 4,563 human-compatible
pairs after the stringent filters. Median counts per library fell from 1,440
to three. The mapping and taxonomic summaries are staged in
[`figures/human_screen/`](../../figures/human_screen/), with no candidate
sequences, read identifiers or human-genome coordinates. The final cohort
report contains all 213 libraries and the completed taxonomic review.
All source, input and output hashes were checked after staging; the report
records the five successful execution IDs and exact source snapshots.
