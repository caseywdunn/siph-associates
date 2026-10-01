# Phase-3 MAG and viral quality decisions

Approved by Casey Dunn on 2026-09-26, after the per-library cohort assembly
was accepted and before catalog construction. The machine-readable thresholds
are in `config/phase3_catalog.json`. Every number below is regenerated from
accepted outputs by `sbatch scripts/phase3_qc_evidence.sbatch`, which writes
the cited tables to `data/results/phase3_cohort/qc_evidence/` (job `11198332`;
minimap2 2.28-r1209, skani 0.3.2).

The assembly stage reports QC values without thresholds. These rules decide
what enters the pooled catalogs. They do not change any per-library output.

## 1. MAG quality

**Decision.** A MAG is retained if CheckM2 completeness is ≥50% and
contamination is <10% (MIMAG medium quality or better). It must also be placed
in Bacteria or Archaea by GTDB-Tk. Bins with completeness ≥90% and contamination
<5% are labelled *near-complete*. Bins recorded as unassessable (CheckM2 found
no DIAMOND annotations; 4 libraries) are not eligible.

**Evidence** (`mag_quality_grid.tsv`; 549 CheckM2-assessed bins):

| Completeness ≥ \ Contamination < | 5% | 10% | 20% | any |
|---|---|---|---|---|
| 90% | 55 | 59 | 61 | 71 |
| 70% | 81 | 90 | 92 | 114 |
| 50% | 104 | **131** | 177 | 204 |
| 30% | 131 | 166 | 214 | 243 |

The 131 retained bins come from 81 libraries: 99 from reference-depleted
libraries (62 libraries) and 32 from fixed-effort libraries (19 libraries).
Most are small, low-GC genomes (median 1.4 Mb, GC 0.34), consistent with
host-associated symbionts.

**Reasoning.**
- The MIMAG medium threshold is the community standard, so the catalog is
  comparable with other studies and the choice needs no special defense.
- Allowing 10–20% contamination would add 46 bins. At that level a bin is more
  likely a strain mixture or chimera, and those errors would propagate into
  competitive mapping.
- CheckM2 models only prokaryotes. Reference-free assemblies retain host and
  other eukaryotic sequence, so GTDB-Tk placement in a prokaryotic domain guards
  against eukaryotic bins receiving prokaryotic quality estimates.
- *Near-complete* is used instead of MIMAG *high quality*, because that tier
  also requires rRNA and tRNA genes, which have not been assessed per bin.

## 2. MAG dereplication

**Decision.** Cluster retained MAGs at ≥95% ANI with ≥50% alignment fraction of
both genomes (the GTDB species convention). The representative is the member
with the highest completeness − 5 × contamination.

**Evidence** (`mag_species_clusters.tsv`, `mag_ani_bands.tsv`; skani over the
131 retained bins):

| ANI / alignment fraction | 95 / 10 | 95 / 30 | **95 / 50** | 97 / 50 | 99 / 50 |
|---|---|---|---|---|---|
| Clusters | 54 | 56 | **59** | 61 | 68 |

Of 930 reported pairs, only 54 fall between 90% and 97% ANI (33 at 90–95, 21 at
95–97). The rest are below 90% (238) or above 97% (638).

**Reasoning.**
- The ANI distribution has a clear gap spanning the conventional 95%
  species boundary, so the cluster count is insensitive to the exact cutoff:
  54–61 clusters across 95–97% ANI.
- The alignment-fraction requirement matters. At 10%, two large Church2025
  clusters (27 and 20 bins) merge through partial alignment. Requiring 50%
  keeps them separate, which is the GTDB practice.
- Completeness − 5 × contamination is the standard representative score (as in
  dRep). It favors complete, clean genomes as mapping references.

**Outcome (2026-09-27).** GTDB-Tk excluded 24 of the 131 quality-selected bins,
all from reference-free libraries, because they carry none of the 120 bacterial
marker genes. They are eukaryotic bins that CheckM2 cannot recognize, and the
evidence estimate above counted them. The catalog therefore holds 107 MAGs in
38 species by greedy centroid clustering (35 by single linkage).

**Amendment (2026-09-30, approved by Casey Dunn): catalog v2.**
- **Rule.** A pair is the same species at ≥95% ANI when *at least one* genome
  of the pair is ≥50% aligned (`align_fraction_rule: at_least_one`), instead of
  both.
- **Why.** The Mycoplasmatales phylogeny showed MAGs at 98.4–99.5% ANI kept as
  separate species. An incomplete MAG cannot cover 50% of its more complete
  partner, so requiring both fractions split a species whenever one member was
  partial. In catalog v1 this left 10 same-species representative pairs
  (7 MAG–MAG), including the 138-library *Physalia* symbiont. Reads were
  split between them, which lowers breadth and inflates breadth-ratio
  (divergent-strain) calls.
- **The earlier concern still holds.** The 10% sensitivity row above merged
  large Church2025 clusters through partial alignment of *both* genomes; the
  amended rule still requires 50% coverage of one genome, and greedy
  clustering is unchanged.
- **Result.** 107 MAGs in 32 species; single linkage also gives 32, so
  representatives no longer link species that greedy clustering separates.
  Six v1 species merge: three clade-A *Physalia* Metamycoplasmataceae into one,
  two pairs in clade B, one DT-91 pair, and one *Cognatishimia* pair.
- **Versioning.** The v1 catalog is kept in `phase3_catalog/mags/v1/`; all
  downstream results are rebuilt as catalog v2 with v1 kept under `v1/`.

## 3. Viral inclusion

**Decision.** A geNomad viral contig enters the vOTU catalog only if it has at
least one geNomad viral hallmark gene **and** is ≥5 kb long or has CheckV
completeness ≥50%. Completeness supported only by terminal repeats (CheckV DTR)
never qualifies a contig by itself. Contigs aligning over ≥50% of their length
to a host reference (or the *N. septata* congener reference for reference-free
*Nanomia*) are reported as **endogenous candidates**, separately from associate
viruses. So are Preplasmiviricota (polinton/adintovirus-like) and Duplornaviricota
(dsRNA virus) contigs without independent evidence of a free virus.

**Evidence: CheckV tiers are unreliable here** (`virus_checkv_tiers.tsv`).
In fixed-effort (reference-free) libraries:
- All 355 CheckV "Complete" contigs rely on DTR completeness. They have a
  median length of 4,962 bp, **none** carries a viral hallmark gene, and 353
  lack any geNomad taxonomy.
- The 376 "High-quality" and 1,243 "Medium-quality" contigs are likewise
  mostly unclassified (361 and 1,148) and rarely carry hallmarks (11 and 84).

**Evidence: host-genome test** (`virus_host_alignment.tsv`). Viral contigs
were aligned with minimap2 (asm20) to the *N. septata* genome:

| Query set | Subset | Contigs | ≥50% aligned to *N. septata* |
|---|---|---|---|
| Reference-free *Nanomia* libraries | CheckV Complete | 50 | 44 (88%) |
| | CheckV ≥ medium | 210 | 120 (57%) |
| | all | 3,491 | 1,087 (31%) |
| | ≥1 hallmark | 905 | 152 (17%) |
| Reference-free non-*Nanomia* | all | 13,908 | 450 (3%) |
| Reference-depleted *N. septata* | all | 44 | 0 |

The reference-free *Nanomia* libraries are congeners of the reference, so the
test can detect host origin there. Most CheckV "Complete" contigs are host
sequence. The controls behave as expected: 0 of 44 contigs from host-depleted
*N. septata* libraries align, and only 3% of distant-host contigs align. For the
other reference-free hosts, the same host signal must exist but cannot be
detected for lack of a reference. That is why the rule relies on positive viral
evidence rather than host subtraction.

**Evidence: rule sensitivity** (`virus_rule_sensitivity.tsv`):

| Rule | Contigs | Libraries | Reference-free / reference-depleted | Unclassified |
|---|---|---|---|---|
| geNomad default | 20,812 | 180 | 17,399 / 3,413 | 16,333 |
| CheckV ≥ medium | 2,342 | 145 | 1,974 / 368 | 2,109 |
| ≥1 hallmark | 3,891 | 177 | 2,276 / 1,615 | 108 |
| **≥1 hallmark and (≥5 kb or ≥50% complete)** | **656** | **134** | 292 / 364 | 1 |
| ≥1 hallmark and (≥10 kb or CheckV ≥ medium) | 370 | 115 | 154 / 216 | 1 |

The selected set is dominated by Uroviricota (322; tailed phages),
Preplasmiviricota (185), Duplornaviricota (69), and Cressdnaviricota (60)
(`virus_selected_phyla.tsv`).

**Reasoning.**
- Host repeats and endogenous elements lack viral hallmark genes but can
  assemble as terminal-repeat "circles". Requiring a hallmark gene removes the
  artifact class identified above and applies identically to reference-bearing
  and reference-free hosts.
- The ≥5 kb or ≥50%-complete floor follows common MIUViG-based practice for
  vOTU inclusion. It excludes gene fragments too short to cluster or map
  reliably. The stricter ≥10 kb variant would roughly halve the catalog without
  addressing a distinct failure mode.
- Even hallmark-bearing contigs can be integrated elements: 17% of those from
  reference-free *Nanomia* align to the congener genome. The plan (D0.4)
  requires host-reference EVEs to be reported separately. Polintons and
  adintoviruses are frequently endogenous in animals, and RNA-virus-like
  sequences in genomic DNA libraries suggest integration. These are therefore
  flagged unless independent evidence supports a free virus.

## 4. vOTU clustering

**Decision.** Cluster included viral contigs at ≥95% average nucleotide
identity over ≥85% of the shorter sequence (MIUViG), using CheckV's
`anicalc`/`aniclust` procedure, after applying the inclusion rule.

**Reasoning.** This is the accepted field standard for species-level vOTUs, so
it needs no data-driven tuning. Clustering after inclusion keeps artifacts from
recruiting genuine viruses into shared clusters.

## Resolved: the recurrent 5,000-bp *Physalia* element is a host retrotransposon

A single near-identical element was assembled at exactly 5,000 bp in 75
Church2025 *Physalia* libraries (21-mer Jaccard 0.93–1.0;
`physalia_5kb_element.tsv`). MEGAHIT marks it circular at about 310× coverage,
and its 141-bp terminal repeat is the k141 circular overlap. It does not align
to the *P. physalis* reference, which is why it survived host depletion. It has
no geNomad taxonomy or hallmark gene, but CheckV called it "Complete" from the
terminal repeat.

Identification (2026-09-26; files in `qc_evidence/`):
- **Nucleotide.** NCBI megablast against core_nt found no significant
  similarity (`physalia_5kb.blastn.ncbi.tsv`), so the element is novel at the
  nucleotide level.
- **Gene content.** prodigal-gv predicts one 1,519-aa ORF covering 91% of the
  ~4.86 kb circle (`physalia_5kb_orf2.faa`). NCBI CD-Search places reverse
  transcriptase (aa 582–760), RNase H (856–989), and integrase (H2C2 zinc
  finger 1091–1146; rve core 1160–1249) domains in that order
  (`physalia_5kb_orf2.cdd.tsv`).
- **Protein similarity.** NCBI blastx against nr: all top hits are Gag-Pol
  polyproteins encoded in animal genomes, at about 39% identity and E = 0
  (molluscs, corals including *Paramuricea clavata* "Transposon Ty3-I Gag-Pol
  polyprotein" and *Pocillopora*, brachiopods, crustaceans, hagfish). There
  were no viral hits (`physalia_5kb.blastx.ncbi.txt`).

**Conclusion.** The domain order RT–RNase H–integrase in a single polyprotein,
with nearest relatives in animal genomes, identifies a Ty3/Gypsy-type LTR
retrotransposon of the *Physalia* genome. It is present at high copy number but
missing from the reference assembly, whose repeat content is incomplete. It is
host sequence, not an associate. The viral inclusion rule already excludes it
(no hallmark gene). The finding supports the decision that terminal-repeat
completeness alone never qualifies a viral contig, and shows that reference
depletion does not remove host repeats absent from the reference.
