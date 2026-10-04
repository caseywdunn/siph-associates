# Mycoplasmatales gene content: results

Analysis per `docs/mycoplasmatales_gene_content_plan.md`, including its
2026-10-02 deviation (relaxed KO tier, calibrated on curated RefSeq genomes).
Tables: `data/results/mycoplasmatales_gene_content/summary/`. Relaxed
assignments are primary; strict-threshold tables (`*_strict.tsv`) are the
sensitivity view. Comparison sets are the near-complete (≥90% complete, <5%
contamination) GTDB r220 Metamycoplasmataceae representatives (n = 130) and
near-complete host-associated external genomes (n = 7).

## Genomes

| Lineage (host) | Genome | Complete % | Contam. % | Size (Mb) | GC | Genes |
|---|---|---|---|---|---|---|
| Clade A (*Physalia*) | MAGSP0005 | 96.2 | 0.0 | 0.82 | 0.297 | 804 |
| Clade B (*Physalia*) | MAGSP0007 | 97.4 | 0.4 | 0.97 | 0.274 | 896 |
| Clade B (*Physalia*) | MAGSP0010 | 95.8 | 0.5 | 0.85 | 0.275 | 780 |
| Clade B (*Physalia*) | MAGSP0029 | 66.8 | 2.3 | 0.57 | 0.278 | 684 |
| Clade B (*Physalia*) | MAGSP0031 | 87.6 | 8.9 | 0.77 | 0.348 | 779 |
| DT-68 (*Nanomia*, *Resomia*) | MAGSP0011 | 91.9 | 0.0 | 0.54 | 0.269 | 552 |
| *Mycoplasma_K* (*Nanomia*) | MAGSP0012 | 91.9 | 0.0 | 0.48 | 0.266 | 443 |

All are small, AT-rich genomes typical of the family (near-complete GTDB
mean 0.87 Mb, 713 genes). Every genome uses genetic code 4: 70–78% of genes
carry an in-frame UGA (Trp), as in *Mycoplasmopsis bovis* (77%), and no contig
with ≥10 genes lacks one (`scripts/gene_content_genome_checks.py`).

`MAGSP0031` has higher GC (0.35) and 8.9% CheckM2 contamination. Its contigs
are uniform in GC (0.31–0.37) and coding density (0.91 overall), and no contig
of ≥10 genes lacks in-frame UGA. The genome is therefore code-4 Mollicutes
throughout, and the contamination estimate more likely reflects duplicated
markers than foreign DNA. Under the plan it still supports presence only.

## Energy metabolism

- **Core glycolysis** (M00002) is complete in every lineage except
  *Mycoplasma_K* (40%; probably missing from the 92%-complete MAG). It is
  complete in 97% of the family, so it is not distinctive.
- **Glucokinase** (K25026) completes the glucose-to-pyruvate module (M00001)
  in both *Physalia* clades. It is found in 25% of family genomes, many of
  which phosphorylate glucose on uptake by PTS instead. These calls are
  relaxed-tier, scoring about 0.9 of the KO threshold. Under strict
  thresholds no genome has M00001 complete.
- **Pyruvate to acetate with ATP.** Both *Physalia* clades encode the 2-oxoacid
  dehydrogenase E1 (scored as the branched-chain paralogs K00166/K00167, see
  plan deviation), E3 (K00382), phosphotransacetylase (K00625) and acetate
  kinase (K00925). This is the PDH–Pta–AckA route found in most family members
  (Pta–AckA module complete in 63%).
- **DT-68** (*Nanomia*/*Resomia*, `MAGSP0011`, 91.9% complete) has no
  detected PDH, Pta or AckA and carries one **L-lactate dehydrogenase**
  (K00016) assignment in both tiers, consistent with lactate fermentation.
  It is the only siphonophore MAG with a strict-threshold K00016 assignment.
  `MAGSP0031` (clade B) also has one K00016 assignment in the primary relaxed
  tier, but none in the strict tier. Because that genome has 8.9% estimated
  contamination, this is a threshold-dependent, single-genome presence
  observation, not an established clade-B trait. Thus lactate dehydrogenase
  is not exclusive to DT-68 under the primary assignment rule.
- **Arginine deiminase pathway** (arcA K01478, arcB K00611, arcC K00926;
  identical in both tiers):
  - present in clade A, DT-68 (three carbamate kinase copies) and
    *Mycoplasma_K*;
  - **not detected in clade B**, in either near-complete genome, although 40–45%
    of family genomes carry it.
  Clade B is therefore expected to depend on sugar fermentation for ATP.
- **Glycerol.** Glycerol kinase (K00864) is in both *Physalia* clades; clade A
  and `MAGSP0031` also have the uptake facilitator (glpF, K02440; both tiers).
  No siphonophore genome encodes
  glycerol-3-phosphate oxidase (glpO, K00105) or dehydrogenase (K00111). The
  H2O2-producing glycerol pathway, a virulence factor in some *Mycoplasma*,
  is absent, as it is from all 130 family comparison genomes.

## Host-derived nutrients

- **Sialic acid.** Only clade A encodes N-acylglucosamine-6-phosphate
  2-epimerase (nanE, K01788; both tiers) and N-acetylneuraminate lyase (nanA,
  K01639; relaxed tier only). These are found in 21% and 13% of family
  genomes. No sialidase (K01186) was detected, so free sialic acid would have
  to come from the host or other microbes.
- **CoA from pantothenate** (M00120) is complete in clades A and B and in DT-68
  (relaxed tier). It is complete in 12% of the family, the distinguishing step
  being the bifunctional coaBC (K13038). Under strict thresholds it is 50–67%
  complete everywhere.
- No chitinase or chitobiase, and no urease.

## Defence and genome maintenance

- No CRISPR-Cas genes in any siphonophore genome (Cas1/Cas9 in 30% of the
  family).
- Restriction-modification: type I (hsdRMS) only in `MAGSP0010` and type III
  (mod/res) only in `MAGSP0007`; none in clade A, DT-68 or *Mycoplasma_K*,
  against 68–87% of family genomes. R-M and CRISPR loci often sit in repeats
  that MAG assembly misses, so these absences are tentative.
- No fatty-acid synthesis (fabH, accA), as expected for Mollicutes;
  cardiolipin synthase is present throughout.

## Single-genome features of `MAGSP0031` (clade B)

Thiamine salvage (thiD/thiE/thiM on contig `k141_808`), complete C1-unit
interconversion (formate–THF ligase, K01938), and sulfur reduction (K18367) are
found only in `MAGSP0031`, on code-4 contigs typical of that genome. They meet
the plan's "rare in family" rule (≤5% of family genomes) at lineage level, but
because they come from one genome flagged for contamination they are reported
as features of that genome, not of clade B.

## Interpretation

For the seven siphonophore MAGs, relaxed and strict tiers agree on the
arginine deiminase, glycerol and nanE assignments. They also agree on the
lactate dehydrogenase assignment in DT-68, but the additional assignment in
`MAGSP0031` is relaxed-only. Concordance here refers to the siphonophore
MAGs; comparison-family frequencies can depend strongly on the tier (for
example, K00016 occurs in 68.5% under relaxed assignments and 10.8% under
strict assignments). The siphonophore Mollicutes share
the reduced, fermentative metabolism of the family and do not encode the
H2O2-producing glycerol pathway, urease or sialidase associated with
pathogenic relatives. The *Physalia* clades ferment sugar to acetate. Clade A
can also catabolize arginine and sialic acid, while clade B appears to lack the
arginine deiminase pathway and so may rely on host-derived sugars. The
*Nanomia* lineages differ: DT-68 is a lactate fermenter with arginine
catabolism.

These are predictions from conserved orthologs. Surface proteins, adhesins
and lineage-specific genes are mostly unannotated by KOfam (plan rule 5).
