# Preparing sequences for public deposition

The deposit workflow makes curated copies of the MAGs and packages the viral,
non-host SSU and trematode-marker sequences. Viral associates and sequences
reclassified as endogenous candidates are packaged separately using
the current reporting manifests under `figures/manuscript_revision/catalog/`.
The package must match their disjoint ID sets and sequence checksums exactly. It preserves the catalog used for
read mapping. The manifest distinguishes candidate GenBank genomes from genomes
retained for an archive. Mixed bins and unreviewed sequence structures are held
separately for provenance only, with unresolved metadata and screening requirements.

The MAG decisions are in `config/mag_deposit_decisions.tsv`,
`config/mag_contig_removals.tsv`, `config/mag_adaptor_trims.tsv` and
`config/mag_sequence_repairs.tsv`. They contain the reviewed evidence paths and
checksums. Runtime rules read these versioned tables and the catalog FASTAs;
they do not depend on the temporary review directory.

`config/mag_organism_names.tsv` maps each MAG's GTDB classification to a proposed
NCBI organism name and a verified taxid for screening. Exact NCBI name lookups
are recorded separately. A missing generic name is marked for NCBI to create;
it is not assigned the taxid of a different organism. The email preparation
table is `config/ncbi_name_verification_request.tsv`. No email is sent by the
workflow.

## Prepare and inspect the sequence bundle

Run from the analysis repository root on McCleary:

```bash
sbatch -p ycga scripts/deposit_preparation.sbatch --dry-run
sbatch -p ycga scripts/deposit_preparation.sbatch
```

The bounded workflow is `workflow/deposit_preparation.smk`. It applies the
reviewed sequence edits, runs CheckM2 on edited MAGs, packages the sequences and
screens all five sequence sets with FCS-adaptor. Job scratch and reports remain
on GPFS. The launcher records the code revision, working diff and source archive.

Inspect `data/results/deposit/bundle/`:

| File | Contents |
|---|---|
| `mag_manifest.tsv` | Every MAG, sequence hash, quality, routing and blockers |
| `curated_mag_catalog.tsv` | Updated MAG quality and size, with original catalog identity |
| `genbank_candidates.tsv` | Genomes meeting the sequence/quality selection rules, pending final screens |
| `archive_genomes.tsv` | Below-threshold genomes retained for the sequence archive |
| `held_genomes.tsv` | Mixed bins or unresolved palindromes retained solely for provenance, outside submission candidates and the claim-supporting archive |
| `mimag_biosample_metadata.tsv` | One draft metadata row per MAG, with source reads where available |
| `sequence_id_mapping.tsv` | Original identifiers and short unique submission identifiers |
| `mags/` and `sets/` | Curated MAG FASTAs and combined sequence sets for screening |
| `source_evidence/` | Existing source/feature evidence for viral, SSU and marker sequences |
| `validation_inputs/` | Draft MAG `.fsa`/`.src` inputs and explicit tool, template and feature-review gaps |
| `status.json` | Current preparation scope and the outstanding author/validation steps |

The MAG selection uses completeness ≥90%, total length ≥100 kb and contigs
≥200 bp, with explicit exclusions for the unrefined mixed-strain MAGSP0028 and
mixed host/bacterial WAM-Z97902 bin.3. Unreviewed exact palindromes hold a genome
out of both the candidate and claim-supporting archive sets. Held genomes remain
in the all-MAG screening input so their diagnostic reports are retained; they
receive no draft submission source files. The criteria follow the [NCBI genome
FAQ](https://www.ncbi.nlm.nih.gov/genbank/wgsfaq/) (checked 2026-10-24).
An eligible genome is still a candidate until the final screens and author
metadata checks have passed. An absent SRA accession remains blank and flagged.

MAGSP0031 is independently reproducible with
`sbatch -p ycga scripts/revision_magsp0031.sbatch`. Its exact reverse-complement
palindrome is halved only after the original length and sequence symmetry pass.
The remaining catalog sequence and all mapped breadth estimates stay unchanged.

## Run FCS-GX after adaptor review

Inspect all FCS-adaptor reports. Resolve any newly found adaptor contamination
in the decision tables and rerun preparation before enabling the FCS-GX target:

```bash
sbatch -p ycga scripts/deposit_fcsgx_preflight.sbatch
# After preflight succeeds and the reports are reviewed:
sbatch scripts/deposit_fcsgx.sbatch
```

The launcher requests 16 CPUs and 720 GB on `bigmem`; the workflow schedules one
650 GB screen at a time. It reads the existing GPFS database and writes to GPFS.
It does not copy the database to `/tmp` or `/dev/shm`. FCS-GX is excluded from the
ordinary preparation target and also requires `enable_fcsgx=true`. The separate target makes the large-memory resource
request explicit. [FCS-GX instructions](https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart).

After the reports are reviewed, resume validation on the final sequences:

1. Record accepted FCS-GX cleaning decisions and repeat affected screens.
2. Finish the NCBI organism-name crosswalk and confirm it with NCBI before
   registering one MIMAG BioSample per genome.
3. Run `table2asn` with validation and discrepancy reports on each submission
   set, using checked source annotations and feature tables. COI boundaries and
   translation frames require the marker evidence; no feature is inferred merely
   from a FASTA header.
4. Run PGAP taxcheck on MAG candidates and ribovore on the submitted rRNA genes.
5. Populate accession fields after the author's manual submission. Prepare the
   archive genomes for the later Zenodo release.

No file produced before these steps is labeled a validated submission. Final
table2asn/PGAP/ribovore jobs follow FCS-GX report review. The sequence bundle contains no invented accessions, author names,
BioProjects or BioSamples.

A real author submission template (`.sbt`) is required for table2asn. The workflow
does not invent contact details or a submission citation. MAG source tables use
the recorded metadata; viral, SSU and marker source and feature tables still need
sequence-specific review. Draft input files do not establish successful validation.

The NCBI name `Metamycoplasmataceae bacterium` (taxid 2979199) is verified. The
installed GX database, built 2023-01-24, assigns its BLAST division only as
`bacteria`. For those genomes, `config/fcs_gx_division_overrides.tsv` explicitly
sets `prok:mycoplasmas` from their reviewed Mycoplasmatales placement. Their
verified NCBI taxid is retained. This screening parameter is distinguished from
NCBI taxonomy and does not change deposited organism names. All other verified
screening taxids resolve in the installed division table without overrides.

Sequence holds are documented in `config/mag_sequence_holds.tsv`. A hold is an
open sequence-review decision, not a claim that an exact palindrome is itself
an NCBI exclusion criterion. MAGSP0020 meets the numerical thresholds but its
7,234-bp full-contig palindrome remains unedited pending review. Existing mean
read depth does not resolve its structure. The manifest reports numerical
thresholds separately from routing and keeps held genomes as provenance only.
