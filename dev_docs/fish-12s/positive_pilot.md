# Fish 12S positive-control pilot

Started 2026-10-17 on branch `fish-12s-screen`, following Casey's identification
of YPM-IZ-110906 and YPM-IZ-110907 as positives from an earlier regex screen.
The [analysis guide](../../docs/fish_12s_screen.md#positive-control-assessment)
defines the methods and reproduction route. Existing sequence analyses are not
dependencies of these standalone jobs.

## Inputs and selection

`config/fish_12s_known_positive.tsv` was joined from the existing
`figures/cohort/sample_metadata.tsv` and `data/metadata/raw_files.tsv`, using
`library_id`; all five requested/previously positive libraries matched, and all
nine paired lane paths exist. The three additional libraries have retained fish
SSU evidence in `figures/cohort/eukaryote_evidence.tsv`. No metadata was edited.

The mitochondrial RefSeq download is
`https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz`.
The listing dates it 2026-09-03. Response headers are in
`data/results/fish_12s/reference/refseq_mito_2026-09-03/download.headers.txt`;
input and output SHA-256 checksums are in the reference-preparation directory.
SortMeRNA 7.1.0 was installed from conda-forge in an isolated environment; its
exact package URLs are in `sortmerna_environment.txt`. No existing environment
was modified.

## Execution

| Job | Scope | Outputs under `data/results/fish_12s/` |
| --- | --- | --- |
| 13499282_0 | Nine synthetic paired controls, regex/BBDuk; passed | `known_positive_controls/2026-10-17/` |
| 13499291_0–8 | All nine selected lanes, regex and BBDuk | `known_positive/2026-10-17/` |
| 13499309 | Bowtie2 index and competitive mitochondrial BLAST database; completed | `positive_references/2026-10-17/` |
| 13499315_0 | Bowtie2/SortMeRNA synthetic controls; passed | `positive_aligner_controls/2026-10-17/` |
| 13499321_0–8 | BLAST candidate verification, depends on the corresponding 13499291 task and 13499309 | `known_positive/2026-10-17/*/*/alignment/` |
| 13499330_2,5 | Bowtie2/SortMeRNA direct raw-read searches, one full lane per user-identified specimen | `positive_aligners/2026-10-17/` |

Statuses of real-data jobs must be read from their execution records; submission
is not completion. The array concurrency limit for regex/BBDuk and BLAST is four.
No full-cohort screen has been submitted.

Both filter implementations passed the same synthetic controls: fish at the
original motif on either mate/orientation, fish outside the motif, and a
cartilaginous-fish variant were recovered. The frog also matched and must be
resolved downstream. Random reads, motif confined to qualities, poly-G, and ATTT
repeat controls were rejected. BBDuk and SortMeRNA retained the expected mates.
The initial synthetic run has 9 input pairs, 3 regex-positive pairs, and 5 pairs
retained by each whole-gene method. The Bowtie2/SortMeRNA ID checks are saved in
`positive_aligner_controls/2026-10-17/synthetic/known_answer/checks.json`.

Every Python-driven run records its source snapshot, manifest/reference hashes,
raw input size and modification time, commands, runtime, exit status, and output
hashes in `run.json`. Large raw FASTQs are identified by the existing lane
inventory and their recorded sizes/times, not rehashed. Shell stages retain the
launch script, traced commands, tool versions, input/output hashes, and a
completion marker written only on success. Cluster logs are in
`logs/fish_12s/`.

Interpretation remains provisional: these libraries were selected for previous
fish evidence, so they cannot estimate cohort-wide detection rates or specificity.
Passing these controls also does not replace reference-held-out sensitivity
tests. A whole-gene candidate is not itself a fish identification.
