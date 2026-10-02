# Mycoplasmatales gene-content plan

Pre-specified on 2026-10-02, before any annotation was run. Settings are in
`config/mycoplasmatales_gene_content.json`; the stage is
`workflow/rules/mycoplasmatales_gene_content.smk` (target
`mycoplasmatales_gene_content`).

## Question

What do the siphonophore Metamycoplasmataceae encode, and how does this
differ from their host-associated and free-living relatives? The aim is to
describe their likely biology (energy metabolism, nutrient uptake,
biosynthetic capacity, defence) and so constrain whether they look like
commensals, nutritional partners, or pathogens.

## Genomes

- **This study (catalog v2):** the seven Mycoplasmatales species
  representatives, grouped by the phylogeny into four lineages:
  - siphonophore clade A: `MAGSP0005`;
  - siphonophore clade B: `MAGSP0007`, `MAGSP0010`, `MAGSP0029`, `MAGSP0031`;
  - *Nanomia*/*Resomia* DT-68: `MAGSP0011`;
  - *Nanomia* *Mycoplasma_K*: `MAGSP0012`.
- **Host-associated relatives:** the 29 external genomes included in the
  phylogeny (octocoral, coral, jellyfish, crustacean and isopod hosts).
- **Family context:** all 143 GTDB r220 species representatives of
  Metamycoplasmataceae, minus any already present as an external genome.

Every genome is re-scored with CheckM2 so all quality values come from one
tool. Genomes enter the comparison at ≥50% completeness and <10%
contamination.

## Methods

- **Gene calls:** pyrodigal 3.6 with genetic code 4 (UGA = Trp), used by all
  Mycoplasmatales.
- **Function:** KofamScan 1.3.0 against the KOfam release current on
  2026-10-02 (SHA-256 recorded). A protein is assigned a KO only when its score
  meets that KO's adaptive threshold. Each protein takes its best-scoring KO.
- **Pathways:** KEGG module completeness from kegg-pathways-completeness
  1.3.0. A module is "complete" at 100%.

## Rules for interpretation

1. **Absence.** A function is called *not detected* in a lineage only when
   it is missing from every near-complete member (≥90% complete, <5%
   contamination, MIMAG high quality). A lineage with one near-complete
   genome is reported as a single-genome observation with its completeness.
   `MAGSP0029` (67%) and `MAGSP0031` (8.5% contamination) cannot support
   absence statements.
2. **Presence.** One threshold-passing KO in any member is a presence.
   Lineage-level presence of a module uses the best member.
3. **Distinctiveness.** A module is distinctive of a siphonophore lineage when
   it is complete there but complete in <25% of near-complete GTDB family
   representatives, or missing from all near-complete lineage members but
   complete in ≥75% of them. These are descriptive contrasts, not
   statistical tests, because genomes are phylogenetically dependent.
4. **Focal functions.** The pre-specified KOs in the config are reported
   for every genome whatever the module results. They cover arginine
   deiminase, glycerol/H2O2, urease, fermentation end products, pyruvate
   dehydrogenase, sialic acid, chitin, CRISPR-Cas, restriction-modification,
   and lipid synthesis. KO identities are checked against the KOfam
   definitions in the output.
5. **Limits.** KOfam annotates only conserved orthologs. Surface proteins,
   adhesins, gliding machinery and lineage-specific genes are mostly
   unannotated, so their absence from KO tables is not evidence of absence.

## Outputs

`data/results/mycoplasmatales_gene_content/summary/`: `genome_stats.tsv`,
`ko_matrix.tsv`, `module_completeness.tsv`, `lineage_modules.tsv`,
`focal_functions.tsv`.
