# Analysis plan (PLoS One)

This is the analysis plan for the siphonophore-associates paper: a **clean,
symmetric** pipeline that treats all three source studies (Church et al. 2025,
Ahuja et al. 2024, Ahuja et al. 2026) as first-class. It is built fresh in this
repo, keyed off [`manifest.csv`](manifest.csv).

The design was piloted in a private exploratory repo (`20260513_siph_symbionts`);
that pilot **justifies the parameter choices locked here** but is never cited in
the manuscript. The pilot grew incrementally with Physalia as the base and the
Ahuja datasets as additions; this plan removes that path dependence (see §7).

Some parameter values remain provisional pending the last pilot runs — see §9.

For the actionable, restart-safe implementation sequence and goal completion
criteria, see [`EXECUTION_PLAN.md`](EXECUTION_PLAN.md). It also records the
manifest and host-depletion issues that must be resolved before cohort work.

**Phase 0 completed 2026-08-29:** 205 unique libraries, 313 explicit paired
FASTQ inputs, 208 study-provenance records, header-derived sequencing batches,
locked resource/audit inventories, and a verified checksum freeze. No cohort
analysis has started; next is the Phase-1 workflow skeleton.

---

## 1. Principles

1. **Three studies are first-class and symmetric.** The unit of analysis is the
   sequencing library/specimen; the source study (Church 2025, Ahuja 2024, Ahuja
   2026) is a provenance *covariate*, not a hierarchy.
2. **One pipeline, applied identically to every library**, regardless of study or
   whether a conspecific reference exists.
3. **The backbone is trim → screen, applied identically to all libraries; host exclusion
   is downstream, not a front-end.** The k-mer screens run on full trimmed reads (they
   ignore host). Host depletion feeds only assembly: Kraken2 `--classified-out` for most
   species, reference mapping for the two with a genome — the latter doubling as a
   *supplementary* host-subtracted analysis, with a robustness check that the two routes
   agree (§5). This keeps the reference asymmetry out of the primary path and Methods.
4. **The reference catalog is built once, from evidence pooled across all libraries
   of all three studies, then every library is mapped against it** — Physalia against
   the *Nanomia*-derived references and vice versa.
5. **No incremental scaffolding in Methods.** Parameters are locked and justified;
   engineering (retries, path fixes, format bugs) never appears.
6. **All four target domains are first-class**, each with detection → reference →
   validation/quantification, applied uniformly across all libraries:
   **(a) bacteria/archaea, (b) viruses/phages, (c) eukaryotic parasites
   (e.g. trematodes), (d) prey** (metazoan/other eukaryotes captured at collection).
   The shared challenge for (c)/(d) is distinguishing non-host eukaryotes from the
   siphonophore host — **Opisthokonta ≠ host** (the Metazoa/eukaryote signal holds
   prey and parasites, not just host), and Kraken2's *Homo* calls are ~99 % false
   and must be validated by mapping.
7. **k-mer detection is nomination only; mapping breadth is the arbiter.** No associate is
   asserted from a k-mer/containment hit. Screens (sylph/Kraken2/phyloFlash) only *nominate*
   candidates; each is then validated by fetching (or de novo assembling) its reference and
   competitively mapping — a claim stands only at sufficient **`covered_fraction` (genome
   breadth)**, paired with `read_count`/`%identity`. Breadth cleanly separated real
   presence (*Alteromonas* 73 %) from k-mer artifacts (*Endozoicomonas* 0.03 %) in the pilot.
   Guardrails: `verify_taxon` on every fetch (the *Pararhizobium* bird-genome incident);
   validate *Homo* by GRCh38 mapping; contamination is a hypothesis tested from breadth +
   cross-library/batch pattern, never from organism identity alone.

---

## 2. Unified dataset

| Study | Libraries | Taxa | Conspecific chromosome-scale reference |
|-------|-----------|------|----------------------------------------|
| Church et al. 2025 | 151 *Physalia* | 5 *Physalia* spp. (utriculus, physalis, megalista, minuta, sp1) | *P. physalis* (this study) |
| Ahuja et al. 2024 | 32-species genome skim | 32 siphonophore spp. across the phylogeny | none (except the *Nanomia*/*Physalia* specimens below) |
| Ahuja et al. 2026 | ~24 *Nanomia* population | 4 *Nanomia* spp. | *N. septata* (this study) |

- **One pooled manifest**, specimen-level **deduplicated**. NA19 and CWD16 are
  distinct Ahuja 2024 libraries reused in the Ahuja 2026 analysis; NA19 was also
  delivered in both local data trees. Church YPM-IZ-104465 and Ahuja NA22 are
  the same raw library and are represented once. Cross-study memberships live
  in a separate provenance table; each raw library occurs once analytically.
- Metadata per library: study, species, ocean/geography, collection date,
  host-reference-available (P. physalis / N. septata / none), sequencing depth,
  library batch (for the contamination test in §5).
- Host-reference availability is an objective per-species attribute, not a per-study
  choice.

---

## 3. Analytical DAG

**The backbone is trim → screen, applied identically to every library.** Host exclusion is
NOT a front-end step — the k-mer screens don't need it (they ignore host as unclassified),
and applying reference-based removal only to the two ref species would break the symmetry.
Host depletion is a **downstream branch** that feeds only the assembly track (and the
supplementary reference-based analysis for the two species that have a genome).

```
S0 Manifest ─► S1 Trim ─► S2 Screening — ALL 205 libraries, on trimmed reads ─────────┐
              cap 400M→fastp    Kraken2+Bracken, sylph (GTDB+OceanDNA),                │
              (trimmed = temp)  phyloFlash 16S/18S      (no host exclusion here)       │
                    │                                                                   │
                    └─► S3 Host-depleted reads (branch → assembly only):               │
                          • ref species (P. physalis, N. septata): map trimmed → drop   │
                            chromosome/mt hitters → non-host reads                       │
                          • all others: Kraken2 `--classified-out` (the SAME run as the  │
                            S2 Kraken2) → non-host reads                                 │
                        └─► S4 Assembly: MEGAHIT on non-host reads (per library)         │
        ┌───────────────────────────────────────────────────────────────────────────────┤
        ├─► TB Bacteria/archaea : sylph/16S screen → catalog → competitive map → MAGs (MetaBAT2/CheckM2/GTDB-Tk)
        ├─► TV Viruses/phages   : geNomad on S4 contigs → CheckV → derep → host-linkage → map-quantify
        ├─► TE Euk parasites    : 18S/28S (phyloFlash+PR2) + barrnap on S4 contigs + COI/mitogenome → non-host euk → parasite ID
        ├─► TP Prey             : COI/18S barcoding (BOLD/SILVA) + prey mitogenome recovery → zooplankton/other prey ID
        └─► S9 Integration: catalog of associates across host phylogeny (discovery; presence + cross-domain links)
```

- **S0 Manifest** — pooled, deduped, metadata (§2).
- **S1 Trim (universal front end)** — cap to 400 M read pairs (§4), then fastp
  (adapter/quality), on **raw reads for all three studies** (Physalia re-processed from
  raw, not prior BAMs). The trimmed reads are a **large regenerable intermediate written to
  scratch** (`/vast/palmer/scratch/dunn/cwd7/...`, 7.5 TB, 60-day purge) — **not** `temp()`,
  because they are consumed twice: by S2 (screens, early) and again by TB competitive mapping
  *after* the catalog is built (§ storage note below). Only the small non-host reads +
  products persist on `/work`.
- **S2 Screening (all 205, identically)** — read-level: Kraken2/Bracken, sylph (GTDB +
  OceanDNA); SSU rRNA: phyloFlash vs SILVA (**16S + 18S**). Run on the full trimmed reads —
  no host exclusion (safer for discovery: no reference-induced loss of host-similar symbiont
  reads). Feeds every domain and the Track-F catalog screen.
- **S3 Host depletion (downstream branch, feeds assembly)** — **[locked]** reference species
  map trimmed reads to the **entire** conspecific genome (all scaffolds + mito) and drop
  **all** host-mapping pairs (keep both-unmapped); all other species reuse the S2 Kraken2
  `--classified-out` non-host reads (no redundant Kraken2). The reference route is also the
  **supplementary** host-subtracted analysis for the two ref species.
  - **References (downloaded from NCBI):** all *Physalia* → *P. physalis* **`GCA_041430235.2`**
    (Church 2025); all *Nanomia septata* → **`GCA_048301705.1`** (Ahuja 2026).
  - **Map to the whole genome, not chromosomes-only** — unplaced/alt-haplotype scaffolds are
    host, so whole-genome mapping removes them correctly (chromosomes-only kept host
    alt-haplotypes as false "non-host" in Phase 1). Safe because S3a confirms the reference
    carries no co-assembled bacterial/symbiont sequence that would delete real symbiont reads.
  - **S3a Host-reference scaffold audit (once per reference).** Classify the unplaced
    scaffolds of the *P. physalis* and *N. septata* assemblies (geNomad + BLAST) before using
    them for subtraction: (i) confirm no co-assembled bacterial/symbiont sequence sits in the
    reference (which would wrongly remove real symbiont reads), and (ii) report co-assembled
    non-host / host EVEs as findings. Pilot: *N. septata* (`GCA_048301705.1`) — 458 unplaced
    scaffolds **audited clean** of bacteria, carrying only host Polinton/adintovirus EVEs.
    The equivalent *P. physalis* audit is also complete: no bacterial/plasmid
    hallmark or conjugation signal; virus-like calls are host EVE candidates.
    Both audited references are accepted for whole-genome host subtraction.
- **S4 Assembly** — MEGAHIT on the S3 non-host reads of each bacteria/eukaryote-rich library
  (objective inclusion threshold). One assembly per library feeds **all four** domain tracks.

**Domain tracks** (each: detect → reference → validate/quantify, uniform across libraries):

- **TB — Bacteria / archaea.** Catalog built once from pooled S2 evidence (≥2 lines
  agree; sylph ANI≥95 gate; contamination screen). **References come from NCBI fetch OR de
  novo assembly:** named candidates are fetched (`verify_taxon` on every fetch); candidates
  with no clean NCBI reference — the common case for the strong marine-MAG signal (e.g. the
  OceanDNA-like hits at 24–52 % containment) — are **recovered from the S4 assemblies**
  (MetaBAT2 → CheckM2 → GTDB-Tk r220). All members (fetched + assembled MAGs, dereplicated)
  form one catalog; **all** libraries are competitively mapped against it
  (`covered_fraction`+`read_count`+%id; **breadth is the arbiter, not k-mers**).
- **TV — Viruses / phages.** geNomad on the S4 contigs → CheckV (completeness/quality;
  retain ≥ medium-quality) → dereplicate viral OTUs across all libraries →
  **host-linkage** (iPHoP and/or CRISPR-spacer matching to the TB bacterial MAGs) →
  map-quantify across libraries. (Pilot: geNomad recovered 3,020 viral contigs → 112 ≥MQ.)
- **TE — Eukaryotic parasites (e.g. trematodes).** Recover eukaryotic SSU (18S/28S) via
  phyloFlash + PR2 reclassification and barrnap on S4 contigs; extract COI/mitogenome
  markers; **separate non-host eukaryotes from host** (Opisthokonta ≠ host — do not
  collapse Metazoa to host); classify parasites (trematode/other) to the finest level
  attainable; quantify by targeted mapping only where a <5 %-divergent reference exists.
- **TP — Prey.** Same eukaryotic marker recovery (COI/18S; BOLD + SILVA/PR2), aimed at
  metazoan/other prey captured at collection (copepods and other zooplankton); recover
  prey mitogenomes where depth allows; distinguish prey from parasites biologically
  (tentacle-captured prey vs tissue-associated parasites) and from host.
- **S9 Integration** — the primary output is a **catalog of detectable associates** across
  the siphonophore phylogeny: which bacteria/viruses/parasites/prey occur in which hosts,
  each presence breadth-validated (not a quantitative abundance comparison). Then the
  qualitative patterns — host-species / host-phylogeny / geographic associations and
  **cross-domain links** (phage↔bacterial-host co-occurrence; parasite/prey↔host) — read as
  presence/association, with `read_pairs` (the exact per-library depth) reported as a
  detection-sensitivity covariate, not used to equalize effort.

TE and TP share machinery (both are non-host eukaryote recovery + marker ID) and are
distinguished at interpretation; they are listed separately because the biological
question (parasitism vs predation) differs.

### Intermediate files, storage, and rerun behavior (Snakemake)

The trimmed reads are needed **twice** — by S2 (screens) and again by TB competitive
mapping *after* the catalog is built from those screens — so they cannot be single-use
`temp()`. They are large (~5 TB capped) and regenerable, so:

- **Write trimmed reads to scratch** (`/vast/palmer/scratch/dunn/cwd7/...`; 7.5 TB free,
  60-day purge) as **normal outputs, not `temp()`**.
- **Downstream rules consume them via `ancient(...)`** and the workflow runs with
  **`--rerun-triggers mtime`**. This yields exactly the desired semantics:
  present → used; purged & a downstream product is pending → regenerated; purged & all
  downstream products done → **not** regenerated (the DAG only builds what's needed);
  regenerating them never cascades reruns of completed work (`ancient()` ignores their mtime).
- **Persistent products live on `/work`** (~300 GB): non-host reads, small filtered
  competitive-mapping BAMs (CoverM inputs), assemblies, MAGs, catalog, reports.
- Consequence: a **scratch purge never triggers a blanket rerun** — only regenerates reads
  when a genuinely pending job needs them. (Tradeoff of `--rerun-triggers mtime`: editing a
  rule's code won't auto-rerun it; force with `snakemake -R <rule>`. Never target the
  trimmed-read paths directly, or they'll be forced to regenerate.)

Host subtraction is **streamed** (`bwa mem | samtools fastq -f 12`), so no full-genome BAM
is ever written; only the non-host reads land on `/work` (§5). Competitive-mapping BAMs are
filtered to catalog-mapping reads only (small) and can be `temp()` after CoverM.

---

## 4. Locked parameters (justified by the pilot)

| Stage | Parameter | Value | Justification (from exploration) |
|-------|-----------|-------|----------------------------------|
| S1 | fastp | `--detect_adapter_for_pe --dont_eval_duplication` | dup-eval OOMs on deep skims; trimmed output identical |
| S2 primary | Kraken2 host depletion | standard DB, `--classified-out`, eager load (never `--memory-mapping` on GPFS) | host unclassified for non-model siphonophores; classified = non-host |
| S2 suppl. | reference routing | proper-pair ≥ 80 % to conspecific reference | clean bimodal split observed (N. septata 84–95 % vs congeners 31–74 %) |
| S3 | reference genomes (NCBI) | Physalia → *P. physalis* `GCA_041430235.2`; *N. septata* → `GCA_048301705.1` | published chromosome-scale assemblies (Church 2025 / Ahuja 2026) |
| S3 | mapping target | **whole genome** (all scaffolds + mito), drop all host-mapping pairs | unplaced/alt-haplotype scaffolds are host; audited (S3a) free of co-assembled bacteria |
| S3 | compute-ceiling cap | **400 M read pairs** (full depth below; trim only the ~5 % of libraries above) | discovery, not quantitative — cap only stops the giant outliers (up to 1.0 B pairs) dominating runtime |
| S3 | sylph DBs | GTDB-r220 c200 + OceanDNA c200 | — |
| S3 | phyloFlash | SILVA 138.1 NR99 (local build), 16S + 18S | resolved species-level parasite/prey in the pilot |
| S4 | assembly | MEGAHIT, `--min-contig-len 1000` | recovered 95.6 % of the *Alteromonas* genome de novo in NA33 |
| TB | catalog inclusion | ≥2 lines agree; sylph **ANI ≥ 95** mappability gate; taxon-verified fetch | ANI<95 won't recruit reads; verify guard prevents accession-swap errors |
| TB | mapping filter | `bwa mem -M -T60` → `view -q30 -F0x900 -f0x0002` → markdup; CoverM identity≥95, breadth≥75, proper-pairs | suppresses repeat-driven multi-mapping; breadth is the arbiter |
| TB | presence call | `covered_fraction` threshold (paired with `read_count`) | e.g. *Alteromonas* 0.73 = real vs. <0.11 trace/artifact |
| TB | MAG QC/taxonomy | CheckM2 + GTDB-Tk r220 | — |
| TV | virus ID / QC | geNomad on contigs → CheckV; retain ≥ medium-quality | pilot: 3,020 viral contigs → 112 ≥MQ |
| TV | phage host-linkage | iPHoP + CRISPR-spacer match to TB MAGs | links phages to their bacterial hosts |
| TE/TP | euk classification | phyloFlash 18S + PR2 v5; barrnap SSU on contigs; COI vs BOLD | Opisthokonta ≠ host; validate *Homo* Kraken2 FP by GRCh38 mapping |
| TE/TP | targeted marker mapping | only where a < 5 %-divergent reference exists | distant refs drop reads under BWA divergence cutoff (pilot Track D2) |

(Values marked ≈ / "e.g." are provisional pending completion of exploration — §9.)

---

## 5. Host handling, and robustness

- **Screens (all libraries):** run on full trimmed reads, no host exclusion.
- **Assembly input (all libraries):** host-depleted — Kraken2 `--classified-out` for most
  species; reference mapping for *P. physalis* + *N. septata*.
- **Supplementary:** the reference-mapping route for the two ref species is also reported
  as a host-subtracted analysis.
- **Robustness checks that de-path-depend the host handling:**
  1. For the two reference species, show symbiont composition from the reference-free
     vs reference-based routes is congruent.
  2. Genome vs transcriptome reference congruence (Ahuja 2026 already established this
     for *Nanomia*; Church 2025 for *Physalia*) — cite as support for reference choice
     not driving results.

---

## 6. Symmetric catalog construction (the key fix)

Build **one** reference set from pooled cross-study evidence, then map everyone against
it. Concretely this means the published comparison "is the *Physalia* community found in
other siphonophores, and vice versa?" is answered from a single mapping pass against a
single catalog — not from a Physalia catalog with taxon-specific additions. The pilot's
finding (Physalia *Vibrio* community ≈absent in *Nanomia*; *Nanomia* carries a distinct
*Alteromonas*-led community) must be reproduced under this symmetric construction.

---

## 7. Path-dependencies this design removes

- Physalia reused pre-made BAMs (different upstream) → **re-process all studies from raw
  identically**.
- Catalog = "30 Physalia genomes + 10 additions" → **one pooled, symmetric catalog**.
- "Phase 1 / Phase 2", "Track A–F" Physalia-centric framing → **stage-based symmetric DAG**.
- Depth cap applied only to some tracks/samples → **full depth everywhere; a single
  400 M-pair compute ceiling trims only the giant outliers** (not rarefaction).
- Reference-based host subtraction as the backbone → **reference-free primary; reference-
  based supplementary with robustness check**.
- Engineering artifacts (CRLF, OOM retries, SLURM path fixes, empty-BAM skips) → absent
  from Methods by construction.

---

## 8. Recompute required (honesty about cost)

Symmetry is not free:
- **Physalia (151 libs): re-trim + re-deplete from raw reads** — the largest recompute.
- **Rebuild the catalog** symmetrically from pooled evidence.
- **Re-map all libraries** (all three studies) against the single catalog.
- Reusable as-is only where byte-identical to the clean pipeline (e.g. the Ahuja fastp
  trims, the existing assemblies).

---

## 9. Open items to resolve before locking this design

These depend on exploration still in progress and must be settled first:
- **Bacteria (TB): final catalog membership** — pending GTDB-Tk of the RF assemblies and
  the reverse Physalia-vs-new-candidates check; resolve the OceanDNA MAG signal (highest
  sylph containment) by de novo assembly vs reference fetch.
- **Viruses (TV): not yet piloted in Phase-2** — run geNomad/CheckV on the new assemblies
  to confirm the track transfers beyond Physalia; decide host-linkage tool (iPHoP vs
  CRISPR-spacer) and the dereplication/OTU criterion.
- **Eukaryotic parasites & prey (TE/TP): under-explored so far** — the Phase-2 pilot has
  been bacteria-focused. Need to run 18S/PR2 + barrnap + COI on the new libraries to
  confirm parasite/prey recovery across taxa, and settle the parasite-vs-prey and
  non-host-vs-host euk decision rules. Was the strongest gap in the exploration.
- **Assembly inclusion threshold (S4)** — set objectively once non-host yields across all
  studies are known.
- ~~Rarefaction depth~~ **Resolved:** compute-ceiling cap set at **400 M read pairs** from
  the exact per-library `read_pairs` (SRA spots + counted); trims only ~5 % (the giant
  outliers), all else full depth.
- **Contamination test design (TB)** — sequencing batch is now normalized per
  FASTQ pair from the Illumina header's instrument, run, flowcell, and lane
  fields. Test the human/skin/soil candidates (*C. acnes*, *Lawsonella*,
  *Bradyrhizobium* — all low-breadth in the pilot) against cross-library/batch pattern; no
  blank controls exist for these skims, so batch structure + breadth is the criterion.
- **Host-handling robustness result (S2)** — confirm the two routes agree before demoting
  reference-based to supplementary.

---

## 10. Repo layout for the paper pipeline (when built)

```
paper/
  DESIGN.md            this spec
  config.yaml          locked parameters (§4), single source of truth
  manifest.tsv         pooled, deduped library manifest (§2)
  scripts/             one module per stage S0–S8, parameterized, study-agnostic
  results/             stage outputs (untracked)
  METHODS.md           generated from DESIGN.md + config.yaml
```

The exploratory tree is frozen as the pilot and is not modified by the paper pipeline.
