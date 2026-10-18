# COI sensitivity and NA35 audit

This follow-up uses preserved pilot results on `fish-12s-screen`. It does not
scan the real libraries again or change the accepted manuscript analyses.
The [execution guide](../../docs/crustacean_coi_pilot.md#reference-divergence-validation-and-candidate-audit)
and [shared search rules](../../workflow/rules/crustacean_coi_search.smk) specify the method.

## Held-out reference test

Validation job **13500312** completed in 3:17, exit 0, peak RSS 4,165,216 KiB.
It selected one annotated RefSeq COI sequence from each of 23 genera, spanning
all seven target groups, and generated five paired 150-nt fragments per sequence.
All panel records mentioning those genera, including original names, and identical
full sequences were removed: 5,725 exclusions, 139,512 retained records.
[Design and selected accessions](validation/design.json) and
[per-pair results](validation/sensitivity_pairs.tsv) retain the exact test.
Genus labels are not a reconciled phylogeny; errors, abundance and detection
probability in real sequencing libraries are not simulated.

| Target group | Positive pairs | BBDuk recruited | DIAMOND recruited | Nucleotide classification: crustacean |
| --- | ---: | ---: | ---: | ---: |
| Amphipoda | 25 | 24 | 25 | 17 |
| Copepoda | 25 | 23 | 25 | 14 |
| Decapoda | 25 | 25 | 25 | 22 |
| Euphausiacea | 10 | 10 | 10 | 6 |
| Lophogastrida | 5 | 5 | 5 | 0 |
| Mysida | 10 | 10 | 10 | 0 |
| Ostracoda | 15 | 15 | 15 | 3 |
| **Total** | **115** | **112** | **115** | **62** |

Here “target” means the seven groups listed above. The current competitor label
`non_target` also includes crustaceans outside those groups, so it must not be
read as a general exclusion of all crustaceans. Broader survey classification
will need explicit coverage of those additional crustacean groups.

The [machine-readable summary](validation/sensitivity_summary.tsv) separates
recruitment from assignment. Of the 115 true target pairs, 28 were ambiguous,
seven unresolved, and 18 incorrectly classified as non-target by the nucleotide
comparison. The five host/fish/insect controls were recruited but classified
non-target; random and poly-G controls were not recruited. All 12 original
exact/synonymous control classifications also still pass after the ranking fix.

This supports retaining translated recruitment for sensitivity. The nucleotide
assignment should not be used as a hard rejection filter for divergent reads.
The three translated-only positive pairs are a measured sensitivity gain, unlike
the original synonymous control, which both methods recovered.

## Candidate ranking and source labels

The original classifier first discarded alignments below 80% identity or 90 nt,
then ranked the survivors. For NA35 pair `c00000016`, this promoted an amphipod
match (80.612%, 98 nt, 92.4 bits) over a stronger decapod-labelled match
(75.333%, 150 nt, 104 bits). The corrected classifier ranks all competitors first
and applies the evidence requirements to the near-best matches. Regression tests
cover this failure and a near-best non-target competitor below the identity cutoff.

Original reports are unchanged. Revised reports in
`data/results/crustacean_coi/validation/2026-10-17/reclassified/` classify the same
92 pairs as 66 unresolved, 23 non-target, two ambiguous and one crustacean.
That remaining provisional crustacean label is NA35 `c00000008`; `c00000016` is
now unresolved. Neither is an accepted prey identification.

Three leading amphipod-labelled barcodes (AY830439, AY830420, AY830436) were
deposited from **gut contents**, as recorded in their GenBank source features.
Their linked paper explicitly examined diets of the source amphipods:
[Blankenship and Yayanos (2005)](https://pubmed.ncbi.nlm.nih.gov/15723681/).
Those organism labels cannot establish the identity of an ingested sequence.
The other suspect barcode, MH425628, is labelled *Pleistacantha kannu* and has
no gut-content source qualifier. Its label remains under review.

The [GenBank annotation audit](validation/genbank_reference_audit.tsv) and
[retrieval hashes](validation/source_record_receipt.json) identify saved source
records. Whole-mitochondrial matches for the four NA35 reads are short/divergent,
with no convincing crustacean assignment. Reviewed-protein matches identify
cytochrome oxidase homologs; the four barcode references have leading bacterial
matches, although these are too divergent to determine their source organisms.
The actual alignments and organism lookup are retained in
[mitochondrial matches](validation/mitochondrial_na35_best.tsv),
[code-5 protein matches](validation/swissprot_code5_na35_best.tsv),
[code-1 protein matches](validation/swissprot_code1_na35_best.tsv), and
[protein source organisms](validation/swissprot_match_organisms.tsv).

## Local bacterial genomes resolve the NA35 pairs

Both provisional NA35 crustacean candidates are **false positives** for
crustacean prey. The [catalog alignments](validation/na35_bacterial_matches.tsv)
and [audited pair labels](validation/na35_audited_pairs.tsv) supersede their
COI-panel labels:

- `c00000016`: both mates match the *Vibrio echinoideorum* reference
  GCF_024347455.1 across all 151 nt, each with two differences (98.675% identity).
  They align in inward-facing orientation on the same contig with a 341-nt span.
- `c00000008`: the mates match the *Alteromonas genovensis* reference
  GCF_010500895.1 at 97.351% over 151 nt and 92.647% over 136 nt. They map to the
  same contig in compatible orientation with a 460-nt span.

These are close bacterial-reference matches, not species-level identifications
of the reads' source bacteria. Their substantially stronger nucleotide evidence
rejects the crustacean interpretation. Neither pair matches the NA35 assembly at
the audit's BLAST threshold. The suspicious MH425628 barcode also aligns to a
*Photobacterium* catalog genome over 650 nt at approximately 86.6% identity;
its crustacean annotation requires curation.

The completed [nr audit](validation/nr_hits.tsv) supports bacterial oxidase
homologs, including cbb3-type cytochrome c oxidase subunit I. Excluding the
barcodes' own deposited proteins, the crab-labelled MH425628 sequence has a
99.6% match across 231 amino acids to *Photobacterium* WP_299017237.1. The three
gut-content barcodes have bacterial matches at 80.3–95.3% across 232–233 amino
acids. Their own deposited proteins also occur in nr; those self matches do not
validate their organism labels.

All four records are flagged in the [reference review](validation/reference_review.tsv)
for quarantine from target assignment. This is a curation decision for the next
panel, not an in-place change to the preserved pilot databases. The exact
bacterial species is not inferred from the short reads or the largest protein
score; close genomic matches and paired-read geometry resolve the NA35 false
positives more directly.

## Execution evidence and next decisions

Protein competition used 103,053 retained target proteins and 17,251 annotated
RefSeq non-target COI proteins. It does not include unannotated cnidarian
barcodes or bacterial cytochrome-oxidase competitors. With a 200-hit limit,
the 80% amino-acid criterion classified 74/115 target pairs as crustacean,
36 ambiguous and five incorrectly non-target. Twenty-three query lists were
still near-best at the reporting limit.

The follow-up using DIAMOND `--top 1` retained 379,242 alignments across 235 reads.
It classified 65/115 target pairs as crustacean, 42 ambiguous and eight incorrectly
non-target at 80% amino-acid identity. The change can expose additional or
stronger competitors because DIAMOND's target-retention option affects search
heuristics as well as output ([DIAMOND documentation](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#output-options)). At 60% and 70% identity it classified 61/115 as
crustacean, with additional ambiguity. All tested settings gave four non-target,
one ambiguous and two unrecruited negative controls; none was called crustacean.
The [limited](validation/protein_summary.json) and
[all-near-best](validation/protein_all_near_best_summary.json) summaries are retained.

Protein competition reduces false non-target labels relative to nucleotide
competition, but short conserved fragments still do not support reliable broad
assignment for every true crustacean. Lowering the identity cutoff does not
resolve that ambiguity. This is a comparison on the same controls, not an
independent estimate of performance after choosing a cutoff.

Execution records are in `data/results/crustacean_coi/executions/{job}/`, including
source snapshots, actual commands, input hashes, versions, exit status and output
hashes. The broad protein checks use existing November 2023 databases, not a claim
of complete current database coverage.

- **13500312:** completed held-out nucleotide benchmark, whole-mitochondrial and
  reviewed-protein candidate checks, and corrected reporting.
- **13500331:** cancelled the overly expensive sensitive nr search before it
  completed; partial output is not evidence of absence.
- **13500343:** focused fast nr search for close homologs of the two NA35 pairs
  and four suspect barcodes; completed in 20:02, exit 0, peak RSS 7,719,120 KiB.
- **13500387:** protein-level competitive assignment on the same held-out controls;
  completed in 2:01, exit 0, peak RSS 622,280 KiB.
- **13500468:** all-near-best protein comparison completed in 1:06, exit 0,
  peak RSS 577,456 KiB.
- **13500629:** bacterial catalog and NA35 assembly checks completed in 0:14,
  exit 0, peak RSS 3,090,276 KiB. The retained `summarize_local_audit.py` regenerates
  the audited labels, checks mate orientation/coordinates and reproduces their hashes.

All submitted analysis jobs have finished or were explicitly cancelled as noted.
Source identities and completion records are summarized in
[completed executions](validation/completed_executions.json). Eleven focused
tests, Ruff lint/format, shell syntax and all three COI rule-graph checks pass.
The final dry run schedules only revised reports; real-library recruitment and
accepted analyses are not scheduled by the validation entry point.

Full-cohort screening has not started. Preserve the union of recruited pairs,
including reads with uncertain or non-target nucleotide labels, for broader
verification. Retain DIAMOND as well as nucleotide recruitment. Add bacterial
oxidase competitors and curate the target barcodes before fixing acceptance
criteria; the four flagged records are an initial finding, not an exhaustive
audit of the panel. Protein comparison alone did not resolve short-fragment
ambiguity. The intended production scope remains both
panels on every library, using all retained host-depleted reads where available.
