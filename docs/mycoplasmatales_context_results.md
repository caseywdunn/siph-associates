# Mycoplasmatales placement and distribution across cnidarians: results

Three follow-on analyses of the siphonophore Mycoplasmatales (catalog v2
species IDs). Gene content is reported separately in
`docs/mycoplasmatales_gene_content_results.md`.

| Analysis | Stage | Config | Outputs |
|---|---|---|---|
| Phylogeny | `mycoplasmatales_phylogeny` | `config/mycoplasmatales_phylogeny.json`, `config/mycoplasmatales_external.tsv` | `data/results/mycoplasmatales_phylogeny/` |
| Cnidarian 16S database | `cnidarian_16s` | `config/cnidarian_16s.json` | `data/results/cnidarian_16s/` |
| *Cassiopea* 16S | `cassiopea_16s` | `config/cassiopea_16s.json`, `config/cassiopea_16s_runs.tsv` | `data/results/cassiopea_16s/` |

## 1. Phylogenetic placement

**Method.** GTDB-Tk 2.4.1 `de_novo_wf` on bac120 markers (FastTree,
WAG+GAMMA), with GTDB r220 Mycoplasmatales and Acholeplasmataceae as the
outgroup. Added genomes: the seven siphonophore species and 36 external
Mollicutes genomes from cnidarian, crustacean and isopod hosts (NCBI BioSample
host search, 2026-09-30, and Keller-Costa et al. 2022). External genomes
needed ≥50% completeness and placement in the order; 36 of 39 were included.
The tree has 483 tips. Support values are FastTree local supports.

**Siphonophore lineages.** All four fall within one 21-tip clade (support
0.999). It holds the seven siphonophore genomes, the four octocoral DT-68
genomes, and the GTDB genera *Mycoplasma_K* (4), *Mycoplasma_O* (2), DT-68 (2;
one is also an octocoral genome added here), RYZV01 and JAIFTM01:

- **Clade A** (`MAGSP0005`, *Physalia*) is sister to *Mycoplasma_O* (*M.
  mahonii* and one uncultured species; support 1.0).
- **Clade B** (`MAGSP0007`, `MAGSP0010`, `MAGSP0029`, `MAGSP0031`) is a clade
  of siphonophore genomes only (support 1.0). Its species are present in 55
  libraries from all five *Physalia* species and one *Nanomia bijuga* library.
  Its sister is the DT-68 group, but with weak support (0.22).
- **DT-68** (`MAGSP0011`, *Nanomia* and *Resomia*) groups with four octocoral
  symbionts, *Leptogorgia sarmentosa* (two), *Eunicella gazella* and *Swiftia
  exserta*, and one further GTDB DT-68 genome. Support 1.0.
- ***Mycoplasma_K*** (`MAGSP0012`, *Nanomia*) falls among the four GTDB
  *Mycoplasma_K* genomes (*M. todarodis*, *M. marinum* and two uncultured
  species; support 1.0).

**Other cnidarian Mollicutes are distant.** The jellyfish genomes (*Aurelia*,
*Rhopilema*) and most coral genomes, including black coral, deep-sea
*Callogorgia* and coral metagenomes, belong to family MT37, mostly
*Spiroplasma_C*. Ten coral genomes have no family assignment and fall outside
the siphonophore clade. Only the octocoral DT-68 genomes are close to a
siphonophore lineage.

## 2. Cnidarian 16S database (McCauley et al. 2023)

**Method.** ASVs from the six regional libraries of the Cnidarian Microbiome
Database were matched over their full length to the 16S genes of the five
siphonophore species with an assembled gene (`MAGSP0005`, `0007`, `0010`,
`0012`, `0031`), at ≥97% identity. Each ASV is assigned to its best-matching
species, with ties listed together.

**Result.** 98 samples from 23 studies carry a matching ASV.

- **Identical (≥99.5%) matches are always rare.** They occur in 7 hosts for
  `MAGSP0005` and 4 for `MAGSP0031`, at no more than 0.6% of a sample's reads.
  `MAGSP0005`-identical ASVs also occur in seawater from 10 samples (≤0.01%).
- **Abundant matches are related lineages, not the same species.** ASVs at
  97.6% identity to `MAGSP0005` dominate *Acropora loripes* (median 92% of
  reads) and *A. sarmentosa* (87%). Both come from a single study (Chan et al.
  2019, V4 region). ASVs at 97.6% to `MAGSP0031` occur in rhizostome jellyfish
  (*Mastigias*, *Cassiopea ornata*) at ≤0.4%.
- **No cnidarian in the database carries a siphonophore species at
  appreciable abundance.** Short amplicons cannot resolve species, so 97–99%
  matches are read as related lineages.

## 3. *Cassiopea* 16S (Florida Keys)

**Data.** BioProjects PRJNA1020388 (96 Illumina V3–V4 libraries: gastrovascular
cavity 36, exterior 36, water 8, sediment 9, swab and filter controls 6, and
one internal sample) and PRJNA1020446 (PacBio HiFi full-length 16S from three
medusae). Downloaded from ENA with MD5 checks.

**Method.** vsearch 2.29.1: pairs merged, primers trimmed (341F/805R), reads
with expected errors >1 removed (860,221 reads kept), denoised with UNOISE3
and chimera-checked (3,703 variants). 85% of filtered reads map to a variant at
97%. Variants (over their full length) and full-length HiFi reads (38,067
reads of 1.3–1.7 kb) were matched to the siphonophore 16S genes from 90%
identity upward.

**Result.**
- No variant and no full-length read matches a siphonophore Mycoplasmatales
  16S gene at ≥90% identity. Florida *Cassiopea* does not carry these
  lineages.
- *Cassiopea* has its own Mollicutes. The closest variants to our genes
  (`ZOTU2`, `ZOTU7`, `ZOTU8`; about 86% to `MAGSP0005`, a distant lineage)
  make up a median 18% of gastrovascular-cavity reads (up to 92%; >1% in 27 of
  36 samples). They are essentially absent from exterior tissue (≤0.1%) and
  were not detected in water, sediment or controls.

## Interpretation

The siphonophore Mycoplasmatales belong to a marine host-associated part of
Metamycoplasmataceae. It also includes octocoral (DT-68) and other
marine-animal (*Mycoplasma_K*, *Mycoplasma_O*) lineages, but not the
*Spiroplasma*-like Mollicutes of corals and true jellyfish. The *Physalia*
clade B is known only from siphonophores. Database and *Cassiopea* surveys
found related lineages in other cnidarians but no other host with a
siphonophore species at meaningful abundance. *Cassiopea* shows the same
pattern with its own lineage: a Mollicutes confined to the gastrovascular
cavity and absent from the environment.

**Limits.** The siphonophore-only status of clade B rests on the genomes and
amplicon data available. The DT-68 sister relationship to clade B is weakly
supported. Amplicon matches below 98.7% over full-length 16S do not establish
species identity.
