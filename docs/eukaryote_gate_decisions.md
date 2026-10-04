# Eukaryotic parasite and prey decisions

Approved by Casey Dunn on 2026-09-30. Rules are frozen in
`config/eukaryote_gate.json`. The evidence is regenerated from accepted outputs:
`workflow/rules/eukaryote_gate.smk` (target `eukaryote_evidence`) produces the
read-lineage table, and `python3 scripts/calibrate_eukaryote_reads.py` writes the
calibration tables to `data/results/eukaryote_gate/calibration/`.

## Evidence sources

- **Assembled SSUs.** Full-length SSUs assembled by phyloFlash (Phase 2), each
  with a SILVA 138.1 best hit giving full taxonomy and identity. 32 libraries
  have none, because phyloFlash's SPAdes step found no rRNA-containing contig.
  This is a sensitivity limit, not a failure.
- **Assembly 18S.** 18S genes found by barrnap on each library's Phase-3
  assembly (429 genes in 165 libraries), classified against the same
  competitive reference.
- **SSU reads.** Each library's phyloFlash-extracted SSU reads, re-mapped
  competitively (below).

## 1. Read-level evidence requires competitive re-mapping

**What went wrong first.** The phyloFlash read mappings failed their negative
control. After correcting a lineage-truncation error in the first extraction,
insect SSU appeared in 182 of 205 libraries at ≥97% identity, and sponges,
annelids, and molluscs appeared almost everywhere. 18S has long stretches that
are nearly identical across animals, and SILVA lacks most of these host
species. Host reads from those stretches therefore landed on arbitrary animal
references.

**Decision.**
- Re-map each library's SSU reads to SILVA 138.1 NR99 plus the cohort's 150
  assembled host siphonophore SSUs (minimap2 `-ax sr -N 200 -p 0.99`,
  secondaries kept).
- Keep pairs whose mates both align over ≥100 bp.
- Assign each pair to the lowest common ancestor of both mates' near-tied hits.
  Reads from conserved regions then resolve only to broad ranks, and host reads
  resolve to the host.
- Exclude from read-level evidence:
  - the host;
  - all Cnidaria, because non-siphonophore cnidarian assignments still occur
    in every library (up to 155,111 pairs at ≥99%), since not every host
    species has an assembled SSU;
  - pairs resolved only to Metazoa, its container clades, or its ancestors
    (Amorphea, Obazoa, Opisthokonta, Holozoa, Choanozoa).

**Result.** After re-mapping, the negative-control background fell to
insects 22 libraries (maximum 36 pairs), land plants 33 (maximum 13), and
non-human mammals 27 (maximum 18). No human reads were found.

## 2. Presence grades

The counts in this section describe the original 2026-09-30 calibration
snapshot. Use the corrected gate summary described in section 7 for reporting
the current results; evidence-row counts and nonredundant detections differ.

| Grade | Rule | Detections |
|---|---|---|
| Validated (reads) | ≥37 read pairs at ≥97% identity, one above the maximum negative-control background (36) | 54 |
| Validated (assembled) | A non-host SSU assembled by phyloFlash or on the library's Phase-3 assembly | 49 (phyloFlash) |
| High confidence | Both, for the same reporting unit and library | 42 |
| Trace | 2–36 read pairs at ≥97%; reported, not claimed as presence | 235 in 104 libraries |

Reporting units are phylum for animals and the first five SILVA ranks for other
eukaryotes.

**Evidence** (`threshold_grid.tsv`, `negative_controls.tsv`,
`assembled_positives.tsv`):

| Minimum pairs (≥97%) | Detections | Negative-control detections | Assembled SSUs recovered |
|---|---|---|---|
| 2 | 295 | 40 | 51/52 |
| 10 | 98 | 9 | 49/52 |
| 25 | 60 | 1 | 46/52 |
| 37 | 54 | 0 | — |
| 50 | 49 | 0 | 41/52 |

**Reasoning.**
- Casey Dunn rejected presence from assembled SSUs alone in order to keep
  sensitivity. Read-level evidence counts toward presence, with its threshold
  set just above the background from implausible lineages, the same logic as
  the bacterial decoys.
- Lowering the threshold to 25 would add 6 detections at the cost of one
  negative-control false positive.
- Novel or divergent lineages have little read support at ≥97% identity to
  SILVA, but they enter through the assembled tier. Examples: a copepod SSU
  at 78–85% identity, a dinoflagellate, and *Paramoeba* at about 87%.
- Weaker read signals are kept as trace rather than discarded.

## 3. Naming depth

Identity to the reference sets naming depth, not presence:
- ≥97%: genus;
- 90–97%: family or order;
- below 90%: novel lineage at class or phylum.

SILVA species labels are never reported, because SILVA labels are unreliable
here: one trematode best hit is labelled "*Daphnia galeata*".

## 4. Verification

Assembled non-host SSUs for the key taxa are checked with NCBI BLAST: the
trematodes, *Paramoeba* and its kinetoplastid endosymbiont, the nematodes,
copepods, and the fish.

## 5. Role field and contamination

- Role is recorded as an evidence field (biological interpretation, not a
  classifier bin):
  - **parasite:** obligate parasitic lineages, for example Digenea, parasitic
    nematode orders, Syndiniales, Apicomplexa, and Kinetoplastea;
  - **prey:** free-living planktonic metazoans, for example copepods,
    chaetognaths, larval fish, molluscs, and annelids;
  - **unassigned:** everything else.
- The lineage lists are in the config.
- Human sequence is flagged as contamination. Single-library detections,
  including the one fish, are kept (Casey Dunn); NCBI verification is the
  check.

**Amendment (2026-09-30, approved by Casey Dunn).**
- **What was wrong.** The parasite list originally included all of
  Kinetoplastea. That labelled the 12 kinetoplastid detections as parasites of
  the host, but 11 of them have NCBI top hits to the "PLO of *Paramoeba
  invadens*" (*Perkinsela*), an obligate endosymbiont that lives inside
  *Paramoeba*. It co-occurs with *Paramoeba* in these libraries.
- **Why lineage alone cannot fix it.** SILVA places the PLO within
  Prokinetoplastina/*Ichthyobodo*, so lineage cannot separate it from the
  parasitic kinetoplastids.
- **The amendment.**
  - Kinetoplastea is replaced by Trypanosomatida in the parasite list.
  - Detections whose NCBI top hit is the PLO are labelled *endosymbiont of
    Paramoeba*.
  - The one read-only kinetoplastid detection, which has no assembled
    sequence to verify, is *unassigned*.
- **Effect.** No grade changes.

## 6. Deferred

- COI and mitogenome assembly, which would allow species-level prey
  identification.
- PR2. SILVA plus NCBI verification is sufficient at the SSU level.

## 7. Reporting correction (2026-10-03)

The manuscript review found two reporting errors. The user authorized their
correction and regeneration of affected downstream outputs. No read alignment,
assembled sequence, reference database, presence threshold, or role decision
changes in this correction.

**Overlapping reporting units.** A short LCA lineage can be an ancestor of a
second reported lineage in the same library. In `Church2025__YPM-IZ-110469`,
324 pairs resolve to `Archaeplastida;Rhodophyceae;Florideophycidae` and 961 pairs
to its descendant ending in `Rhodymeniophycidae`. Both exceed the presence
threshold, but the ancestor reads do not establish an additional organism.

The grader therefore preserves every evidence row and its original grade,
read count, assembled support, and interpretation, and adds:

- `count_as_detection`: true only for validated or high-confidence evidence
  without a validated descendant reporting unit in the same library;
- `reporting_status`: `detection`, `unresolved_ancestor`, or `trace_evidence`;
- `overlapping_descendant_units`: the validated descendants explaining why
  an ancestor cannot be counted separately.

Ancestor reads are not assigned to descendants or used to increase their
grades. A trace descendant cannot suppress a validated ancestor. Ambiguous
ancestor evidence remains available in the full supplement table and can be
shown separately in figures. Count only `count_as_detection == "true"` when
summarizing detected reporting units. These counts are not species richness:
reporting resolution differs among clades and one unit can contain multiple
organisms.

**Naming ceiling versus resolved rank.** The previous `naming_depth` field
said `genus` whenever identity was at least 97%, including classifications
ending at the order `Calanoida`. Identity does not restore taxonomic
resolution lost through LCA assignment. The replacement `naming_ceiling`
records only the maximum naming resolution permitted by the identity rule.
`lca_lineage` and `lca_terminal_taxon` retain the actual selected LCA result;
`named_lineage` retains the conservative name used before this correction.
An order-level LCA with a genus naming ceiling remains an order-level result.

**Validation and resulting counts.** The bounded regrade was checked against
all original rows: every pre-existing evidence field except the renamed
`naming_depth` field is unchanged. All 295 rows remain, including 62 validated
evidence rows and 233 trace rows. Removing the one redundant validated ancestor
from detection summaries gives **61 nonredundant reporting-unit detections in
39 libraries**, including 42 high-confidence detections. Their roles are
12 parasite, 9 prey, 11 endosymbiont of *Paramoeba*, and 29 unassigned. Role
labels remain biological interpretations, not evidence of interactions with
the siphonophore host.

`tests/test_eukaryote_reporting.py` checks overlapping ancestors, multiple
descendants, separate libraries, trace descendants, unchanged read counts,
coarse LCA names, and rejection of an incorrect counting flag. The validator
now checks the reporting annotations as well as the original grade rules.
Accepted pre-correction outputs and producing code are preserved with a
SHA-256 manifest at
`data/results/eukaryote_gate/review_archive/2026-10-03_before_reporting_correction/`.

Regenerate only the affected summaries from the analysis repository root:

```bash
python3 -B workflow/scripts/grade_eukaryotes.py \
  --config config/eukaryote_gate.json \
  --reads data/results/eukaryote_gate/read_lineages.tsv \
  --assembled data/results/eukaryote_gate/assembled_nonhost.tsv \
  --verification data/results/eukaryote_gate/ncbi_verification.tsv \
  --output data/results/eukaryote_gate/eukaryote_grades.tsv
python3 -B workflow/scripts/validate_eukaryote_gate.py \
  --config config/eukaryote_gate.json \
  --grades data/results/eukaryote_gate/eukaryote_grades.tsv \
  --assembled data/results/eukaryote_gate/assembled_nonhost.tsv \
  --verification data/results/eukaryote_gate/ncbi_verification.tsv \
  --output data/results/stages/eukaryote_gate.done
```

The upstream competitive mapping and NCBI verification do not need rerunning
for this reporting-only correction.
