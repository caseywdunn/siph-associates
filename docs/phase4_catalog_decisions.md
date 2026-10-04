# Phase-4 catalog decisions

Approved by Casey Dunn on 2026-09-27, before any Phase-4 catalog was built.
Machine-readable settings are in `config/phase4_catalog.json`. Nomination
evidence is regenerated from accepted outputs with
`python3 scripts/phase4_nomination_evidence.py`, which writes the cited tables
to `data/results/phase4_catalog/nomination_evidence/`.

The catalogs are frozen as version `v1` before cohort mapping. Any later
addition creates a new version and requires a complete remap.

## 1. Bacterial and archaeal nomination

**Decision.** Nominate every genome detected by sylph at ≥95% adjusted ANI in
any library. That is 471 genomes: 333 GTDB r220 species representatives and
138 OceanDNA species representatives, from 117 libraries. Bracken (≥10 reads)
and phyloFlash (≥1 read) support at the same genus in the same library is
recorded per genome as evidence for grading, not used as a filter. This
replaces the PLAN.md rule that at least two lines must agree.

**Evidence** (`sylph_nominations.tsv`, `gtdb_corroboration.tsv`,
`gtdb_uncorroborated.tsv`):

| Source | Hits | Genomes | Libraries | Median ANI |
|---|---|---|---|---|
| GTDB r220 | 654 | 333 | 103 | 97.08 |
| OceanDNA | 278 | 138 | 82 | 97.56 |

- Genus-level corroboration covers 559 of 654 GTDB hits (85%): Bracken 543,
  phyloFlash 473.
- The 64 uncorroborated GTDB genomes are mostly uncultured marine lineages whose
  genera exist only as GTDB placeholders or are absent from the NCBI-based
  Kraken database: SW10, *Nitrosopelagicus*, *Luminiphilus*, UBA8752,
  *Pseudocolwellia*, UBA9214, CAJXJP01, UBA12191, …
- Their sylph matches are as strong as the corroborated ones (median ANI 96.75
  versus 97.14).

**Reasoning.**
- A name-matching agreement rule can only confirm taxa that NCBI and SILVA
  name. It would systematically exclude the uncultured marine lineages that
  PLAN.md expects to carry the strongest signal.
- sylph's ≥95% ANI containment match is specific, and it is the mappability
  gate itself: a genome below it would not recruit reads.
- Catalog inclusion is not a presence claim; mapping breadth is the arbiter.
  An extra genome in a competitive catalog is cheap and reduces misassignment
  of reads to relatives.
- Recording agreement as evidence preserves the information for the evidence
  grades.
- Common reagent and skin contaminant genera remain in the catalog (96 hits,
  15 genomes, 52 libraries; `contaminant_genera.tsv`: *Cutibacterium* 31,
  *Pelomonas* 21, *Bradyrhizobium* 18, *Mesorhizobium* 14, …). The
  contamination tests (D0.4) need them present so their breadth and batch
  distribution can be assessed.

## 2. References, screening, and dereplication

**Decision.**
- **Sources.** GTDB genomes come from the local r220 package, which holds all
  113,104 species-representative genomes, with checksums recorded. Each GTDB
  accession must be a current, unsuppressed NCBI assembly. OceanDNA genomes are
  extracted from the official species-representatives archive (figshare
  10.6084/m9.figshare.c.5564844, file 29164842, MD5-verified). This is the same
  set the sylph OceanDNA database was built from.
- **Screening.** Every reference passes CheckM2 at the MAG thresholds (≥50%
  complete, <10% contamination). OceanDNA genomes must also be placed in
  Bacteria or Archaea by GTDB-Tk.
- **Pooling.** Screened references and the 38 Phase-3 MAG species
  representatives are clustered at ≥95% ANI and ≥50% alignment fraction. The
  representative is the member with the highest completeness − 5 ×
  contamination, whether MAG or reference.

**Evidence** (`mag_overlap.tsv`). All 8 MAG species with GTDB species names
(*Serratia nevei*, *Vibrio caribbeanicus*, *Cobetia amphilecti*, three
*Pseudoalteromonas*, *Vibrio alginolyticus*, *Photobacterium angustum*) are
also nominated by sylph, so MAGs and references overlap. The other 30 MAG
species have no reference and enter only as MAGs.

**Reasoning.**
- One quality bar for every genome prevents weaker references entering. The
  OceanDNA references are themselves MAGs and can be contaminated.
- One genome per species avoids splitting reads between near-identical genomes,
  which would lower the breadth of both.
- Choosing by quality rather than origin usually prefers complete isolate
  genomes, which improves breadth estimates, and keeps novel MAGs as their
  species' only genome.

**Amendment (2026-09-28, approved by Casey Dunn before any mapping).**
- **Substitution.** A GTDB accession that is not current in its own NCBI
  database is replaced by a verified current substitute: the identical current
  record in the other database (GenBank for a suppressed RefSeq record), or the
  latest version of the same accession. The substitute is downloaded from NCBI
  with its published MD5, must be ≥99% ANI to the GTDB genome it replaces, and
  passes the same CheckM2 screen. Without a verified substitute, the reference
  is excluded.
- **Why.** Under the first build, 22 references were excluded, including
  well-supported species such as *Vibrio neptunius* and *V. tasmaniensis*
  (3 libraries each) and *Pseudoalteromonas atlantica* (9 libraries, 99.3%
  ANI). Losing them would send their reads to relatives or leave them unmapped.
- **A status-check bug.** 13 of those 22 were current GenBank records that the
  check could not match. NCBI returns the paired RefSeq record for a GenBank
  query, and the first implementation keyed status only by that record's
  accession. The check now judges each accession in its own database.
- **Result.** 383 of 393 GTDB accessions are current. The remaining 10 have
  substitutes: 7 identical GenBank records for suppressed RefSeq entries and
  3 newer versions.

**Amendment (2026-09-30, approved by Casey Dunn): catalog v2.** Pooling uses
the amended Phase-3 species rule: ≥95% ANI with ≥50% alignment of at least one
genome of the pair (see `docs/phase3_qc_decisions.md` §2). Requiring both
fractions left 10 same-species pairs among v1 representatives (7 MAG–MAG),
which split reads between genomes of one species. All of v1 is kept unchanged
under `phase4_catalog/v1/`; catalog v2 requires a complete remap (Phase 5) and
reanalysis (Phase 6), per the versioning rule in EXECUTION_PLAN.md.

**v2 outcome (2026-10-02).** Catalog v2 has 448 non-decoy genomes: 24 MAG
representatives and 424 references. Including 20 decoys gives 468 mapping
targets. The viral catalogs are unchanged. The distinction between 448 real
targets and 468 total targets was clarified during publication review on
2026-10-03; no catalog sequence or mapping changed.
- **A mixed bin.** The v1 genome present in 138 *Physalia* libraries
  (`BAC00033`, a 54%-complete MAG in 486 contigs) combined the clade-A
  Metamycoplasmataceae with unrelated sequence. 110 of its contigs align to
  the complete clade-A genome (`MAGSP0005`) and are 92% protein-coding. The
  other 376 are 20% coding, have no BLASTN (core_nt) or BLASTX (nr) hits for
  the contigs tested, and do not align to the *P. physalis* assembly
  (GCA_041430235.2). In the 94 libraries where v1 called `BAC00033` present but
  v2 calls the clade-A species absent or trace, 2,030,019 of 2,031,231 reads
  (99.9%) on `BAC00033` fell on those 376 contigs (median 0.03% on the
  *Mycoplasma* contigs per library, at most 0.5%). This sequence recruits reads
  across all five *Physalia* species. It is most likely host sequence absent
  from the assembly, which therefore escapes host depletion. In v2 the bin
  merges into the complete genome, and only 17 of the 376 contigs align to any
  v2 genome. The clade-A species is present in 44 libraries at median breadth
  0.85, not 138 at 0.34.
- **Screen of the remaining MAGs.** Contigs <50% coding (Prodigal) make up
  ≤3.7% of any v2 MAG representative (`BAC00116`, *Xenococcus*, 3.7%;
  `BAC00006`, 2.3%; all others about 1% or less). None align to the *Physalia* assembly.
  Because presence needs ≥10% breadth, such contigs cannot produce a presence
  call by themselves. They add about 5% of reads to `BAC00116` and 4% to
  `BAC00006`. No catalog v3 is needed.
- **Reproduce:** `scripts/catalog_mixed_bin_evidence.py` (run under SLURM,
  48 GB; writes `data/results/phase4_catalog/v2/bin_evidence/`). The BLAST
  queries (2026-10-01, RIDs `BY66C9AB014`, `BY6F2EWT016`) were run
  interactively and are not scripted.
- **Low-confidence lineage.** `BAC00028` (67%-complete Metamycoplasmataceae,
  `MAGSP0029`) is present in 9 libraries at median breadth 0.12, with 8
  divergent-strain calls, and is reported as low confidence.

## 3. Presence-rule controls

**Decision.** Build the controls into the frozen catalog. Lock the presence
rule from their mapping behaviour before any presence result is interpreted.
- **Positives.** Each retained MAG must pass in the library it was assembled
  from.
- **Negatives.** 20 decoy genomes from lineages implausible in marine animal
  tissue are included and mapped under identical filters in all 205 libraries.
  They are GTDB r220 species representatives drawn with seed 20260927 from
  Bacteroidaceae (4), Lachnospiraceae (4), Bifidobacteriaceae (2),
  Streptomycetaceae (3), Pasteurellaceae (3), Nostocaceae (2), and
  Methanobacteriaceae (2). Any decoy within 95% ANI of a catalog genome is
  replaced.
- **Starting rule to evaluate** (D0.2): ≥10% genome breadth, ≥95% identity,
  MAPQ ≥30, proper pairs.

**Reasoning.**
- D0.2 and the Phase-4 acceptance require the presence rule to be locked from
  positives, negatives, and decoys. A breadth threshold without controls cannot
  be defended.
- MAGs are guaranteed positives in their source libraries, so a rule that fails
  them is too strict.
- Decoys compete with real genomes under identical filters. Their breadth
  distribution measures the apparent presence produced by conserved-gene
  recruitment and cross-mapping alone.
- Drawing decoys from gut, soil, and mammalian-mucosal lineages across the
  catalog's phyla (including an archaeal family) samples that artifact across
  lineages, because conserved-gene recruitment is lineage-dependent.

## 4. Viruses and host linkage

**Decision.**
- **Catalogs.** The viral mapping catalog is the 319 associate vOTU
  representatives from Phase 3. The 153 endogenous-candidate clusters are a
  separate catalog, used for reporting only.
- **CRISPR host linkage.** Spacers are called with MinCED in every bacterial
  and archaeal catalog genome and matched to vOTU representatives by BLAST. A
  link requires the full spacer length with ≤1 mismatch.
- **iPHoP** is deferred.

**Reasoning.**
- A CRISPR spacer that matches a virus records a past infection. It is the
  most specific host-linkage evidence available and runs locally on genomes
  already in the catalog.
- The full-length, ≤1-mismatch criterion follows common practice for spacer
  matches.
- iPHoP needs a database of several hundred GB and integrates indirect signals.
  It can be added later if CRISPR links are sparse, as they may be for
  Mycoplasmatales and other host-associated lineages with few CRISPR arrays.
- Host linkage does not feed the mapping catalogs, so deferring iPHoP delays
  nothing else.

## 5. Eukaryotic parasites and prey

**Decision.** Handle these in a separate gate after the bacterial and viral
catalogs are frozen, with its own evidence review.

**Reasoning.** The decisions differ in kind:
- marker databases (PR2 version, COI references);
- separating host, prey, and parasite, which is a biological interpretation;
- marker-appropriate breadth rules.

Bundling them would delay catalogs that are ready now. The evidence (phyloFlash
SSU, barrnap markers on 195 assemblies, and the endogenous-candidate set) is
preserved for that gate.
