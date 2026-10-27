# Submitting the recovered sequences to NCBI

This guide connects this project's prepared sequence files to NCBI registration,
validation and submission. NCBI guidance was checked on **2026-10-26**. Download
current portal templates when submitting; the local metadata tables are drafts.
The sequence deposit is not yet submission-ready, and no NCBI submission has
been made by this workflow.

Use [sequence preparation](deposit_preparation.md) for the reproducible cleaning
and packaging steps and [local validation](deposit_validation.md) for tool setup
and checks. Run their commands from this repository's root on McCleary. Those
commands prepare files; they do not register records or upload sequences.
Completed checks need repeating when their sequence inputs change, not merely
because the submission guide has been updated.

## Progress and remaining work

This checklist reflects the local reports inspected on 2026-10-26:

- [x] Prepare the curated bundle: 104 MAGs, initially routed as 56 GenBank
  candidates, 44 archive genomes and four provenance-only holds.
- [x] Screen all five sequence sets with FCS-adaptor; no findings were reported.
- [x] Run FCS-GX on all 104 MAGs; four reports contain actions requiring review.
- [x] Run ribotyper on 65 non-host SSUs and 37 MAG 16S sequences; all passed.
- [x] Run draft table2asn checks on 56 candidates. This is completion of a draft
  check, **not successful final submission validation**.
- [ ] Resolve all applicable FCS-GX findings and sequence holds; regenerate the
  final candidate list and repeat checks for changed sequences.
- [ ] Verify missing source-read accessions and complete collection metadata.
- [ ] Confirm the submission route, BioProject arrangement and organism names
  with NCBI; then register the required records.
- [ ] Complete the planned PGAP taxonomy check and choose the annotation route.
- [ ] Supply real submitter details and finish source/feature review for each
  sequence class; inspect final validation reports.
- [ ] Submit the reviewed files, resolve NCBI feedback and record accessions.
- [ ] Verify released records and update the paper and sequence archive links.

Evidence is under `data/results/deposit/`: `bundle/`, `fcs_adaptor/`,
`fcs_gx/`, `executions/fcsgx-13932897/` and
`executions/sequence-validation-13934591/`. The latter contains ribotyper results
and `table2asn_drafts/status.json`, which records `submission_ready: false` and
the absent submitter template. The preparation-time `status.json` files and
manifest blocker strings do not incorporate all these later checks; use the
execution reports and their input checksums to assess the current files.

The initial 56-candidate table is not an upload list. In particular,
`Ahuja2026__NA31__bin.4` (`SIPH0003`) has an unresolved FCS-GX `REVIEW_RARE`
finding. Six MAGs lack source accessions in the bundle, including candidates
`Church2025__YPM-IZ-111839__bin.1` and
`Church2025__YPM-IZ-111241__bin.1`. Proposed local accession matches still need
verification against primary archive records. PGAP has not been installed or
run through the proposed deposit workflow.

## 1. Confirm the archive route before registering records

These assemblies reuse reads from several earlier studies. NCBI distinguishes
assemblies of submitters' own sequence data from third-party assemblies of
public data. Ask NCBI which route applies to each source study and sequence
class, how to link the existing projects and reads, and what release conditions
apply. The TPA route has publication-dependent release and excludes single
genes, so it cannot be assumed to cover every marker here.
See [NCBI's Third Party Assembly guidance](https://www.ncbi.nlm.nih.gov/genbank/tpa/).

The author handling submission should contact `genomes@ncbi.nlm.nih.gov` for
MAG routing and taxonomy review; use `gb-admin@ncbi.nlm.nih.gov` for unresolved
GenBank marker/viral routing. Include the source-study/SRA crosswalk, the planned
sequence classes, and the organism-name table. Save the response and resulting
decisions. No email is sent by the repository's workflow.

| Prepared material | Local bundle file | Submission route to confirm |
| --- | --- | --- |
| Final eligible bacterial MAGs | `mags/<mag_id>.fna`; candidate `.fsa` drafts under `validation_inputs/genbank_candidate/` | Genome submission, one genome per file and one organism-specific MIMAG BioSample per MAG; confirm primary/third-party handling |
| Viral associates (401 sequences) | `sets/viruses.fna` | GenBank sequence submission, separate from the bacterial genome batch |
| Endogenous-virus candidates (163 sequences) | `sets/viral_endogenous.fna` | Separate review of host-derived/integrated sequence status and appropriate source annotation before choosing a route |
| Non-host SSU rRNA sequences (65) | `sets/ssu.fna` | GenBank rRNA workflow, subject to the reused-read routing decision |
| Trematode markers (12) | `sets/markers.fna` | GenBank marker workflow appropriate to each sequence; use COX1 only for reviewed mitochondrial COI records |
| Below-threshold MAGs and held sequences | `archive_genomes.tsv`, `held_genomes.tsv` | Separate sequence archive; held sequences remain provenance-only, outside the claim-supporting genome set |

These are paths relative to `data/results/deposit/bundle/`. The combined
`sets/mags.fna` is a screening input, not a batch upload file. Viral associates
and endogenous candidates have disjoint reporting manifests; preserve those
sets and their original-to-submission identifier mapping.
The [Genome portal](https://submit.ncbi.nlm.nih.gov/about/genome/) directs viral,
phage and single-locus sequences to
[GenBank sequence submission](https://submit.ncbi.nlm.nih.gov/about/genbank/).

## 2. Finalize the candidate sequences and source evidence

Use `mag_manifest.tsv`, `genbank_candidates.tsv`, `sequence_id_mapping.tsv` and
`mimag_biosample_metadata.tsv` in the bundle alongside the versioned
[deposit decisions](../config/mag_deposit_decisions.tsv),
[sequence holds](../config/mag_sequence_holds.tsv) and
[organism-name crosswalk](../config/mag_organism_names.tsv).
Apply decisions to the curated deposit copies, preserving the original mapping
catalog and analysis identifiers. Record final file hashes and the code revision.

The current NCBI MAG gate is CheckM/CheckM2 completeness of at least 90% and
assembly length of at least 100,000 nt. Each MAG must represent one organism's
reconstructed genome. See the [MAG instructions](https://www.ncbi.nlm.nih.gov/genbank/wgsfaq/#metagen).
The corrected MAGSP0031 has 87.49% completeness and remains below that gate.
Apply eligibility to the final edited sequences; a scientific quality category
does not itself establish submission eligibility.

Review FCS-adaptor and FCS-GX reports, document decisions for every applicable
action, and screen any changed sequences again. A completed screen does not
resolve a mixed bin or an unexplained sequence structure. Keep tool/database
versions, input hashes and reports with the submission evidence. See
[NCBI FCS](https://github.com/ncbi/fcs).

Verify each source run against its archive record and the original library
identifier. Keep source-study, host specimen, tissue, source BioSample and run
accessions distinct from the new bacterial MAG records. Do not infer a match
from host name alone or invent missing collection data. Missing accessions and
conflicting records remain explicit blockers until resolved or an alternative
is agreed with NCBI.

## 3. Register the BioProject and MAG BioSamples

Use the BioProject arrangement confirmed in step 1. Before MAG BioSample
registration, send NCBI the stable isolate IDs and proposed organism names or
unmodified GTDB lineages for verification. Start with
[`config/ncbi_name_verification_request.tsv`](../config/ncbi_name_verification_request.tsv),
filtered to the final submission set. Screening taxids and FCS division overrides
are not substitutes for accepted organism names.

Create one MIMAG BioSample per MAG, including separate MAGs from the same host
specimen. Preserve the stable `SIPH####` isolate IDs in the bundle. The host/read
BioSample belongs in the source crosswalk, not in place of a MAG BioSample.
Follow the [MAG-specific instructions](https://www.ncbi.nlm.nih.gov/genbank/wgsfaq/#metagen).

Download the current [MIMAG batch template](https://www.ncbi.nlm.nih.gov/biosample/docs/packages/MIMAG/)
from [BioSample submission](https://submit.ncbi.nlm.nih.gov/subs/biosample/).
Transfer reviewed values from the draft metadata table into that template;
do not upload the internal table unchanged. Complete the required environmental
attributes, including collection date, geographic location, coordinates,
isolation source and environmental context, using the source records and
permitted missing-value terms. Preserve documented host and tissue information.

For MAGs, set `sample_type=metagenomic assembly`; put the contributing
SRR/ERR/DRR run accessions in `derived_from`, comma-separated. Follow this
MAG-specific convention even where generic BioSample help mentions SAMN.
Also include the source runs in each contig's `[SRA=...]` FASTA modifier.
See the [MAG instructions](https://www.ncbi.nlm.nih.gov/genbank/wgsfaq/#metagen).
Record the returned BioProject and MAG BioSample accessions in a versioned
crosswalk and propagate them into submission metadata.

## 4. Choose genome annotation and validate the final payload

NCBI accepts unannotated FASTA and can run PGAP during processing, or accepts
annotated ASN.1 `.sqn` files. Local PGAP taxonomy checking and NCBI annotation
are different steps. The current project plan still has the former pending.
For this project's submitter annotation route, preregister locus-tag prefixes
through BioProject/BioSample before preparing features. Local PGAP output has
a separate prefix-remapping option during NCBI processing. See the
[Genome portal annotation instructions](https://submit.ncbi.nlm.nih.gov/about/genome/)
and [local PGAP submission instructions](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/).

For FASTA submission, provide author/contact and submission citation information
in the portal. For the local table2asn route, create a real
[submitter template](https://submit.ncbi.nlm.nih.gov/genbank/template/submission/)
and supply the reviewed source and feature files. The existing
`scripts/deposit_sequence_validation.sbatch` runs structural drafts without
that template; it is not a final-package launcher. The underlying
`workflow/scripts/run_table2asn_draft_checks.py` accepts `--submitter-template`
but remains a draft-check helper, not an uploader or proof of readiness.

Inspect table2asn validation and discrepancy reports, not just its exit status.
Resolve Error/Reject validation messages and FATAL discrepancy findings, and
document warning dispositions. The present draft reports
include unresolved taxonomy/source issues. Use the
[table2asn instructions](https://www.ncbi.nlm.nih.gov/genbank/table2asn/)
and [validation explanations](https://www.ncbi.nlm.nih.gov/genbank/genome_validation/).
Repeat sequence-dependent validation after any cleaning or feature change.

## 5. Submit the approved genome batch

Once the earlier steps are complete, use the
[Genome submission portal](https://submit.ncbi.nlm.nih.gov/subs/genome/)
with the route NCBI confirmed. For a standard batch:

1. Choose Batch/multiple and the correct WGS status for the fragmented MAGs.
2. Supply a current Genome Info table linking each exact filename to its MAG
   BioSample, assembly method/version, genome coverage and sequencing technology.
   Use recorded assembly/depth evidence, not a guessed coverage value.
3. Upload one file per genome. Keep contigs separate within it, with unique
   sequence identifiers; do not concatenate unrelated contigs. Genome contigs
   must be at least 200 bp. Assign circularity, chromosomes or plasmids only
   where supported.
4. Group only compatible records: one BioProject, one file format, common
   WGS status, release date and gap policy; maximum 400 genomes per batch.
5. Review the rendered metadata, release choice and validation messages before
   submitting. Save the submitted files, checksums and portal receipt.

See the [genome submission guide](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/).
If NCBI specifies third-party handling, follow that route instead of assuming
the standard form and release options apply unchanged.

## 6. Submit viral and marker sequences separately

After resolving the reused-read routing question, use the applicable
[GenBank sequence workflow](https://submit.ncbi.nlm.nih.gov/help/genbank/).
Review source descriptors and feature coordinates for each record using
`data/results/deposit/bundle/source_evidence/` and the corresponding analysis
outputs. Those evidence files are not finished submission feature tables.

For SSUs, check organism assignments, rRNA features, orientation and partial
ends; ribotyper PASS alone does not settle taxonomy. For trematode markers,
review which sequences contain rRNA or mitochondrial COI, and verify COI
boundaries, reading frame, genetic code, translation and partial ends.
For confirmed trematode mitochondrial COI, NCBI lists
[translation table 21](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG21).
The portal can propose rRNA/COX1 features, which still require review. Separate
loci and submission types; source-table IDs must match FASTA IDs. These loci
were assembled from DNA sequencing, so an rRNA feature does not make the
submitted molecule RNA.

For viral associates, use the Virus workflow and supply supported features;
these novel viral candidates do not fall within the named groups with automatic
annotation. Separate partial and complete sequences and resolve taxonomy with
NCBI where needed. Where a viral BioSample is required, use the current
[MIUVIG host-associated package](https://www.ncbi.nlm.nih.gov/biosample/docs/packages/MIUVIG.host-associated.6.0/)
and document the source, assembly and lack of viral enrichment accurately.
The siphonophore providing the library is not necessarily the organism a virus
infects, especially for bacteriophages. Keep endogenous candidates separate
until their host/integration interpretation and source descriptors have been
reviewed; a viral classifier
match does not establish an independently occurring virus or a complete genome.
See [GenBank sequence submission instructions](https://submit.ncbi.nlm.nih.gov/about/genbank/).

## 7. Record accessions and check release

Maintain a versioned TSV crosswalk from original sequence/MAG IDs through
submission IDs to BioProject, BioSample, assembly and nucleotide accessions.
Include accession versions, submitted-file hashes, dates and release status.
A portal `SUB` identifier is a tracking receipt, not a public sequence accession.
Keep account credentials and private contact details out of this public repo.

Respond to NCBI validation queries and retain the resulting sequence/metadata
changes with their evidence. For ordinary GenBank genomes, release follows the
selected date or publication, whichever comes first; see the
[genome submission guide](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/).
For TPA, confirm the publication requirement and notify NCBI when it is met;
see [TPA release guidance](https://www.ncbi.nlm.nih.gov/genbank/tpa/).

After release, retrieve the records and confirm the sequences, organism names,
source links and accession crosswalk. Then update the manuscript, supplementary
tables and repository data links. Publish the separately reviewed sequence
archive with its manifests and holds clearly identified. The existing v0.5.0
Zenodo DOI is a code snapshot, not the reconstructed sequence deposit.
