# Halistemma reference addition, 2026-10-19

User authorized the mapping pilot, full host depletion and assembly of both
H. rubrum libraries. Scope and commands: [Halistemma analysis](../docs/halistemma.md).
User explicitly requires waiting for both Bargmannia and Halistemma to be
finalized before the single revised-catalog full-cohort remapping.

NCBI reference: GCA_987347065.1, jhHalRubr1.1, chromosome-level, taxid 316183,
released 2026-09-08. NCBI dataset metadata, assembly sequence inventory and
published checksums are preserved under `data/sources/halistemma/`.
The nuclear and mitochondrial sequences are included in the reference.

The original workflow configuration and sample identifiers remain unchanged.
The two-sample entry point uses new results and retained-read roots, and reuses
the existing cohort computational rules without changing the running Bargmannia
implementation. Prior trimming provenance supplies raw-file identity and exact
read-accounting checks. A seeded million-pair pilot checks each population;
poor mapping stops before full depletion. The 50% pair-removal check is an
operational usefulness gate, not a new biological presence criterion.

No cohort-wide remapping was running at initialization. The Bargmannia workflow
ends at per-library results, as does this workflow. The explicit hold is stored
in `config/post_expansion_remapping.json` and linked from both analysis guides.

Validation: four focused tests passed, covering FASTA identity/alphabet,
paired FASTQ IDs and qualities, preserved raw-read checksums/counts, and
pair-versus-read mapping accounting with low-mapping and missing-library stops.
Ruff and launcher syntax checks passed. The production-profile dry run schedules
33 steps for exactly CWD12 and CWD17, with no pooled-catalog or full-cohort
mapping rules. The generated rule graph passes its freshness check.

NCBI's inventory contains 5,466 nuclear sequences (4,378,693,217 bases) plus
the 13,934-base mitochondrial sequence OZ505748.1. The validator accounts for
NCBI's separate primary-assembly and organellar length conventions. Existing
raw-file paths and compressed sizes match the completed fish-screen records;
full checksum agreement is checked after trimming before the pilot proceeds.

Submitted controller **13568531** on 2026-10-19. The initial seven-day controller
request was rejected by Slurm's maximum wall-time policy; the three-day request
was accepted. The launch record is
`data/results/halistemma/2026-10-19/executions/13568531/`, and its controller log
is `logs/halistemma/controller-13568531.log`. The entry point is explicitly
bounded to `halistemma`; no remapping dependency is submitted.

2026-10-20 status: reference verification and both million-pair pilots passed.
Host mapping removed 99.8711% of pilot pairs in CWD12 and 98.5620% in CWD17.
Full CWD12 depletion completed at 04:00, retaining 131,703 of 100,871,052
quality-filtered pairs. This is below the existing 175,000-pair assembly floor;
the assembly eligibility gate will require review once CWD17 finishes. Preserve
the original full-read assembly as a comparison baseline. Do not change the
running depletion job or silently lower the assembly threshold. CWD17 full
depletion remains active under job 13570279. Cohort-wide remapping stays on hold.

## Host-depletion table, 2026-10-20

The [reporting workflow](../docs/host_depletion_report.md) generated
`figures/host_depletion/mapping_stats.tsv` and its Markdown counterpart at 06:11:
163 rows, 162 complete mappings, and CWD17 explicitly pending. Counts use
individual reads after quality filtering, with pair-preserving host removal.
All 161 original-library counts agree with the independent cohort eligibility
table. Eleven focused tests, Ruff, launcher syntax, and graph freshness passed.
The bounded dry run schedules only one reporting rule; a subsequent dry run
does nothing. Requiring all completed inputs correctly fails on the two absent
CWD17 records without scheduling upstream work.

The main analysis Snakefile was left unchanged: its dependency graph can rebuild
sequence analyses, so this cross-workflow report has a dedicated entry point
consuming existing records. No sequence analyses were launched for the report.

Refresh job **13587095** was submitted at 06:13 with
`afterok:13570279`; Slurm confirmed it pending on that dependency. It will require
all 163 completed mapping records and regenerate the same table with CWD17's
final counts. Resources: one CPU, 2 GB, five minutes. Launcher:
`scripts/host_depletion_report.sbatch`; log:
`logs/host-depletion-report-13587095.log`. It does not depend on the subsequent
assembly gate and cannot start full-cohort remapping.
