# Bargmannia expansion, 2026-10-19

Scope and commands: [analysis guide](../docs/bargmannia.md).

The original fish screen completed and was recorded on `fish-12s-screen` at
`d179f24`; the new work is on `bargmannia-expansion`. Crustacean work remains
on its separate paused branch. All original sample and sequence identifiers
remain unchanged.

Public SRA checks now establish release of all eight new runs. The local
compressed-file sizes match the original uploads, including the formalin
T676-SS3 library. Workflow MD5 checks are a prerequisite to processing, so
file identity does not depend on resolving the conflicting voucher numbers.

The source tables agree on the published species names. D725-D8 is
B. marlena; D611-D6 and T411-SS2 are B. pyrosa. Table S1, Table S2, SRA and
local filenames agree on D725-D8; the methods' D725-D9 is retained as an
additional discrepancy note. The earlier unresolved voucher discrepancies
remain recorded without inventing corrections:

- T676-SS3: S1 YPM-IZ-111936 versus S2/SRA YPM-IZ-11936.
- D611-D6: S1 YPM-IZ-112080 versus S2/SRA YPM-IZ-107003.
- Older Ahuja T596-SS12: S2 assigns YPM-IZ-111936, conflicting with S1's
  T676-SS3 record. No old-library voucher has been changed.

Validation before launch: exact eight-run selection; all six published species;
one unique local read pair per accession; all original-upload byte sizes;
recorded tissue for all eight; public BioSample and collecting-event joins.
The checksum gate rejects a changed file even when its size is unchanged.
Both focused tests pass. The workflow dry run contains 150 steps, restricted
to Manko2026 inputs and the new output root. Assembly and screen commands reuse
the established cohort code; the optional assembly eligibility gate is empty
for the original cohort. The production launcher records source, controls,
job logs and completion status. Combined catalog construction and global
mapping follow review of both Bargmannia and the newly authorized Halistemma
host-depleted assemblies; see the explicit hold below.

Launched on 2026-10-19 with controller job **13567924**. All eight
paired-file checks passed against the MD5 checksums of the public NCBI uploads;
this also confirms the local formalin library's identity as T676-SS3. The first
four trimming jobs started successfully. The controller log is
`logs/bargmannia/controller-13567924.log`; the frozen launch source, metadata,
control results and job records are in
`data/results/bargmannia/2026-10-19/executions/13567924/`.

Do not edit the active workflow or its inputs while this controller is running.
The current manuscript still reports the original 205 libraries; the eight new
libraries will enter the manuscript after combined-catalog mapping and review.

2026-10-19 update: user requires **one** revised-catalog full-cohort remapping,
after both Bargmannia and Halistemma are finalized. Do not launch a Bargmannia-only
intermediate remapping. See `config/post_expansion_remapping.json` and
[Halistemma run plan](halistemma-expansion.md). Neither per-library workflow has
an automatic path into cohort-wide remapping.

2026-10-20 status: all eight assemblies, marker/viral recovery, binning and
CheckM2 assessment completed; the assembly validator passed with zero errors.
The assembly summary contains 36 candidate bins, of which three meet the initial
CheckM2 medium-quality screen. These are not accepted MAGs until taxonomic and
domain review. Initial Kraken2/Bracken, sylph and phyloFlash outputs are complete.

The first controller exited after 109/150 steps because the new fish inventory
writer used `r1`/`r2` whereas both shared consumers require `r1_path`/`r2_path`.
All failures occurred before fish sequence screening began. Corrected the
writer, added a consumer-column regression test, and passed all three Bargmannia
input tests and Ruff. The resumed dry run contains 49 fish-only steps (including
eight regenerated inventories and the aggregate target); no trimming, assembly,
initial classifier screens or cohort-wide mapping is scheduled.
Resumed with controller **13586761**, preserving the failed execution record
and recording the corrected source under its new execution ID.
