# Mycoplasmatales placement and distribution across cnidarians: results

The phylogeny places seven siphonophore Metamycoplasmataceae representatives;
16S comparisons ask where related sequences occur in other hosts. The GTDB
*Tyloplasma*-related `MAGSP0033` is outside these comparisons. The original
catalog remains the source for the phylogeny and read mapping; curated deposit
versions and their quality estimates are described in
[deposit preparation](deposit_preparation.md).

| Analysis | Stage | Config | Outputs |
|---|---|---|---|
| Phylogeny | `mycoplasmatales_phylogeny` inputs, then `scripts/revision_iqtree.sbatch` | `config/mycoplasmatales_phylogeny.json`, `config/mycoplasmatales_external.tsv`, `envs/iqtree.yaml` | `data/results/mycoplasmatales_phylogeny/iqtree3/2026-10-25/` |
| Cnidarian 16S database | `cnidarian_16s` | `config/cnidarian_16s.json` | `data/results/cnidarian_16s/` |
| *Cassiopea* 16S | `cassiopea_16s` | `config/cassiopea_16s.json`, `config/cassiopea_16s_runs.tsv` | `data/results/cassiopea_16s/` |

## 1. Phylogenetic placement

**Method.** IQ-TREE 3.0.1 on the GTDB-Tk 2.4.1 bac120 alignment, with
GTDB r220 Mycoplasmatales and Acholeplasmataceae as the outgroup.
ModelFinder selected LG+F+I+R10 by BIC. Five independent searches each used
1,000 SH-aLRT replicates and 1,000 ultrafast bootstrap replicates with NNI
refinement; the best-scoring run supplies the reported tree and support.
See the [IQ-TREE execution guide](mycoplasmatales_iqtree.md).
Candidate genomes: the seven siphonophore species and 32 external
Mollicutes genomes from cnidarian, crustacean and isopod hosts (NCBI BioSample
host search, 2026-09-30, and Keller-Costa et al. 2022). External genomes
needed ≥50% completeness and placement in the order; 29 of 32 were included,
giving 36 added genomes in total.
The tree has 483 tips. Support values below are SH-aLRT / ultrafast bootstrap
percentages.

**Siphonophore lineages.** The seven representatives fall within one 21-tip clade
(support 99.9/100). It holds the seven siphonophore genomes, the four octocoral DT-68
genomes, and the GTDB genera *Mycoplasma_K* (4), *Mycoplasma_O* (2), DT-68 (2;
one is also an octocoral genome added here), RYZV01 and JAIFTM01.

The figure's host and habitat annotations come from
[`config/mycoplasmatales_source_hosts.tsv`](../config/mycoplasmatales_source_hosts.tsv)
and the study specimen metadata. Reference sources were checked against
BioSample records and source publications. They include octocorals,
cephalopods, a basket star, a pteropod, a scorpion and a foraminifer, as well
as sinking marine particles and a hydrothermal microbial mat. Environmental
recovery does not identify an organismal host or establish free-living biology;
these rows have a dash in the Host clade column. The pteropod name
*Limacina rangii* follows the source publication, while its BioSample records
*L. antarctica*; the difference is retained explicitly in the evidence table.
The generated
[`source_phylogeny_hosts.tsv`](../figures/mycoplasmatales/source_phylogeny_hosts.tsv)
links each displayed row to its evidence. Render the retained tree with
`python workflow/scripts/plot_mycoplasmatales_figures.py` on a compute node;
the figure manifest records the source table and tree checksums.

Within this clade:

- **Clade A** (`MAGSP0005`, *Physalia*) is sister to *Mycoplasma_O* (*M.
  mahonii* and one uncultured species; support 100/100).
- **Clade B** comprises `MAGSP0007`, `MAGSP0010` and `MAGSP0029`.
  `MAGSP0031` is a divergent sister lineage, with 91.5% 16S identity and
  about 56% amino-acid identity to the other clade-B genomes. The branch joining
  these genomes has support 100/100; its relationship to DT-68 is weakly supported
  (64.2/71). Related 16S variants occur in other marine hosts.
- **DT-68** (`MAGSP0011`, *Nanomia* and *Resomia*) groups with four octocoral
  associates, *Leptogorgia sarmentosa* (two), *Eunicella gazella* and *Swiftia
  exserta*, and one further GTDB DT-68 genome. Support 100/100.
- ***Mycoplasma_K*** (`MAGSP0012`, *Nanomia*) falls among the four GTDB
  *Mycoplasma_K* genomes (*M. todarodis*, *M. marinum* and two uncultured
  species; support 100/100). It was reconstructed from the single specimen NA31,
  which is excluded from sample-based biological analyses; prey carry-over is
  a possible explanation. The genome is retained for phylogenetic context.

**Other cnidarian Mollicutes are distant.** The jellyfish genomes (*Aurelia*,
*Rhopilema*) and most coral genomes, including black coral, deep-sea
*Callogorgia* and coral metagenomes, belong to family MT37, mostly
*Spiroplasma_C*. Ten coral genomes have no family assignment and fall outside
the siphonophore clade. Only the octocoral DT-68 genomes are close to a
siphonophore lineage.

## 2. Cnidarian 16S database (McCauley et al. 2023)

**Method.** ASVs from the six regional libraries of the Cnidarian Microbiome
Database were matched over their full length to the 16S genes from MAGs in five of seven
siphonophore genome clusters (`MAGSP0005`, `0007`, `0010`,
`0012`, `0031`), at ≥97% identity. Each ASV is assigned to its best-matching
clusters, with ties listed together. The reviewed misbinned 16S copy on
`WAM-Z97916 bin.2` contig `k141_5955` is excluded by
`config/mag_contig_removals.tsv`. This copy belongs to clade B; retaining it
would incorrectly add `MAGSP0005` to tied clade-B assignments.

**Result.** 98 samples from 23 studies carry a matching ASV.

- **Close (≥99.5%) matches are rare.** The
  [sample-level table](../tables/mycoplasmatales/cnidarian_16s_matching_asv_samples.tsv)
  records their host and seawater occurrence. Matches to clade B include
  *Astroides calycularis* and *Aplysina aerophoba*. The 16S genes in
  `MAGSP0007` and `MAGSP0010` are 99.74% identical, so this threshold cannot
  distinguish those genome clusters.
- **Abundant matches are related lineages, not the same species.** ASVs at
  97.6% identity to `MAGSP0005` dominate *Acropora loripes* (median 92% of
  reads) and *A. sarmentosa* (87%). Both come from a single study (Chan et al.
  2019, V4 region). ASVs at 97.6% to `MAGSP0031` occur in rhizostome jellyfish
  (*Mastigias*, *Cassiopea ornata*) at ≤0.4%.
- **The closest matches are uncommon in the sampled communities.** Short
  amplicons cannot resolve species, so these matches identify reference
  affinities rather than confirmed occurrences of the same bacterial species.

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
  16S gene at ≥90% identity in these datasets. This non-detection does not
  establish absence from Florida *Cassiopea*.
- *Cassiopea* has its own Mollicutes. The closest variants to our genes
  (`ZOTU2`, `ZOTU7`, `ZOTU8`; about 86% to `MAGSP0005`, a distant lineage)
  make up a median 18% of gastrovascular-cavity reads (up to 92%; >1% in 27 of
  36 samples). They are essentially absent from exterior tissue (≤0.1%) and
  were not detected in water, sediment or controls.

## Interpretation

The siphonophore Mycoplasmatales belong to a part of Metamycoplasmataceae
recovered from diverse hosts and marine environments. It includes octocoral (DT-68) and other
marine-animal (*Mycoplasma_K*, *Mycoplasma_O*) lineages, but not the
*Spiroplasma*-like Mollicutes of corals and true jellyfish. Database and *Cassiopea* surveys
found related lineages in other cnidarians, while the closest matches were rare.
The sampled Florida *Cassiopea* instead had a more distant Mollicutes lineage,
concentrated in gastrovascular-cavity samples and undetected in the sampled
water, sediment and controls.

**Limits.** The DT-68 sister relationship to clade B is weakly
supported. Amplicon matches below 98.7% over full-length 16S do not establish
species identity.
