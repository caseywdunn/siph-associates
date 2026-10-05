# Direct comparison with Waki's trematode sequences

This bounded follow-up asks whether markers deposited by Waki and colleagues
can be compared directly with existing siphonophore sequence assemblies. The
accepted parasite screen used 18S/SSU rRNA; Waki's deposited data use 28S rRNA,
ITS2-containing amplicons, and mitochondrial COI. Similar reference names in
nonoverlapping markers are not a direct sequence comparison.

## Scope and inputs

Only the eight specimens with accepted trematode SSUs are searched in their
existing Phase-3 assemblies: YPM-IZ-104464, 110436, 110557, 110631, 110694,
110846, 110882, and 110931 (all `Church2025__` libraries). These assemblies
total about 50 MB. The existing accepted helminth SSU collection supplies
12 sequences from 11 libraries, including both the trematode and cestode
sequence in YPM-IZ-110436. It is included as a marker-overlap control and to
look for physical links between an accepted SSU and a newly compared marker.

The original phyloFlash archives were inspected in all eight specimens. They
retain the extracted rRNA sequences, not full SPAdes scaffolds. The installed
phyloFlash implementation extracts SSU intervals using barrnap and removes
the SPAdes directory. Consequently, the full Phase-3 assemblies are the
available source for additional loci. No read mapping or assembly is rerun.

Reference files, accession metadata, feature annotations, and retrieval
provenance are in `data/sources/waki2026/`. Newly deposited study sequences
and older comparison sequences from Waki's supplement remain separate
reference sets. The current panel additionally includes two study records
found in GenBank that were omitted from the article's accession inventory:
LC889189.1 (28S) and LC889237.1 (COI). The final panel contains 108 newly
deposited study records (26 28S, 11 ITS2-containing, and 71 COI) plus 113
external comparison records (43 28S, 22 ITS2-containing, and 48 COI), for
221 sequences. The initial 106-record study panel and
its first analysis are historical exploratory outputs, not the final selected
comparison.

## Selected findings

The expanded search completed on 2026-10-05 in about 12 seconds inside
SLURM job 11793267. It retained 1,980 assembly-reference HSPs. The
906 alignments at least 200 bases long include distant homology and
conserved sequence, and must not be interpreted as 906 parasite detections.
The comparison adds marker affinities within the already selected specimens;
it does not change the accepted library or detection counts.

The clearest comparisons are:

| Specimen | Host | Marker and reference affinity | Exact local alignment |
|---|---|---|---|
| YPM-IZ-110882, Hawai‘i | *Physalia utriculus* | *Dinurus tornatus* ITS2-containing amplicon, LC889198.1 | 651/651 bases identical (100%) |
| YPM-IZ-110882 | *P. utriculus* | *D. tornatus* COI, LC889235.1 | 851/867 (98.15%) |
| YPM-IZ-110931, Japan | provisional *Physalia* sp1 | *D. tornatus* ITS2-containing amplicon, LC889198.1 | 662/665 (99.55%) |
| YPM-IZ-110931 | provisional *Physalia* sp1 | *D. tornatus* COI, LC889234.1 | 851/863 (98.61%) |
| YPM-IZ-110436, Florida | *P. physalis* | *D. barbatus* ITS2-containing amplicon, LC889203.1 | 672/673 (99.85%) |
| YPM-IZ-110631, northeastern US | *P. physalis* | *Prodistomum* Type 4 28S, LC889183.1 | 1,142/1,149 (99.39%) |
| YPM-IZ-110694, northeastern US | *P. physalis* | *Prodistomum* Type 4 28S, LC889183.1 | 1,157/1,167 (99.14%) |

The *D. tornatus* COI references derive from metacercariae in Japanese
*P. utriculus*. The ITS2-containing *Dinurus* references used above derive
from adults in *Coryphaena hippurus*. The Type 4 28S reference derives
from a metacercaria in *P. utriculus*. These are direct comparisons to
homologous deposited sequences, providing more specific hypotheses than
the original 18S-only reference labels.

Alternative references temper those hypotheses. The Hawaiian and Japanese
COI fragments remain much closer to *D. tornatus* than to the newly added
*D. barbatus* COI on the same subject coordinates. Florida has no close
COI match in this assembly; its ITS2-containing sequence is 99.82%
identical to *D. barbatus* versus 99.08% to *D. tornatus* over 542 shared
subject bases. This is a smaller separation than in the Hawaiian/Japanese
COI comparisons. *D. longisinus* lacks ITS2 and COI representation in
this selected panel, so the earlier 18S closest-reference label cannot be
excluded by complete same-marker sampling here.

The two Type 4-like 28S sequences also closely resemble Type 5 references.
Over 945 shared subject bases, Type 4 versus Type 5 identities are
99.26% versus 98.94% in 110631 and 99.05% versus 98.94% in 110694.
Their COI matches to Type 4 are substantially more divergent (about
87–88%), consistent with a broader *Prodistomum*-related affinity rather
than a demonstrated identity to Waki's Type 4.

The remaining selected specimens (104464, 110557, and 110846) yielded
only short or more distant comparisons to this reference panel. Retaining
them explicitly in the 24-row summary prevents a failure to recover a
close additional marker from becoming an undocumented sample exclusion.
All accepted-SSU matches to Phase-3 contigs were only 46–48 bases long,
and the reference-to-SSU control had no alignment longer than 107 bases.
No accepted SSU was physically linked to the additional loci under the
declared criteria. Consequently, report *D. tornatus*-like,
*D. barbatus*-like, and *Prodistomum*-related sequences in specimens,
without assigning an individual parasite or asserting infection location.

## Reproduction

The serial analysis uses the existing `snakemake` conda environment, with
Python 3.11, Biopython 1.81, and BLAST+ 2.5.0. Exact observed versions and
commands are captured by each run. The analysis ran inside an existing
one-CPU, 4-GB SLURM allocation; it does not require a new assembly job.

From the repository root, using a fresh output directory:

```bash
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python \
  scripts/waki_comparison/compare_waki_assemblies.py \
  --include-context --output data/results/waki_comparison/reproduction
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python \
  scripts/waki_comparison/summarize_waki_comparison.py \
  --analysis data/results/waki_comparison/reproduction \
  --output data/results/waki_comparison/reproduction/report
```

For an independent allocation, ensure `logs/slurm/` exists and submit
`sbatch scripts/waki_comparison/batch.sh` from the repository root. Its default
output path contains the job ID. Set `WAKI_RUN_DIR` to choose another fresh
path. The scripts refuse to overwrite previous completed execution records.

## Methods and reporting

`config/waki_comparison.json` declares the existing assemblies, references,
accepted SSUs, and parameters. All references are searched against the pooled
eight assemblies using `blastn -task blastn -word_size 11 -dust yes -evalue
1e-10`, one thread, and up to 1,000,000 target sequences. No per-target HSP
limit is imposed. Separate searches compare the accepted SSUs with those
assemblies and the reference markers with the accepted SSUs.

Every significant HSP is retained with reference and contig coordinates,
alignment length, identity, scores, lengths, and both aligned sequences.
Alignments at least 200 bases long are flagged as long homology candidates.
This flag is neither a presence threshold nor a taxonomic assignment: a
restricted trematode reference panel can also match distant homologues and
conserved regions in other organisms.

The compact report has all eight specimens crossed with the three markers,
including explicit absence of a significant study-reference hit. For each
specimen and marker, the focal match is the highest-scoring local alignment
to a newly deposited Waki reference. Exact matching-base counts supply the
reported percentage; percentages are not calculated by rounding BLAST's
already rounded three-decimal output a second time. The report retains the
reference's host and developmental stage as biological context.

For competing taxa, eligible HSPs must cover at least 80% of the focal
contig interval. All eligible reference haplotypes are then scored on the
same shared, unambiguous subject-base coordinates, before selecting the
best reference of each taxon. This avoids choosing a longer but less
similar haplotype before normalizing the compared interval. Subject
insertions relative to the reference count as nonidentical; reference
insertions have no subject coordinate and are excluded from this separate
comparison. The original gap-inclusive HSP identities remain available.
These normalized identities are diagnostics, not formal species distances.

Waki's paper labels combine isolate-specific GenBank organism names into
the appropriate named species or numbered type. Two older sequence labels
are reconciled only where Waki explicitly provides the identification:
Hemiuridae sp. A with *Dinurus barbatus*, and Sclerodistomidae sp. A with
*Bathycotyle branchialis*. Original accession and organism labels remain
in the detailed output. Different isolates or these historical aliases
must not be counted as independent competing taxa.

ITS2 records annotate the entire 5.8S/ITS2/28S amplicon as one feature,
without internal boundaries. Therefore report their identities as
ITS2-containing amplicon identities. Short matches can be restricted to
conserved flanks and cannot establish a variable-ITS2 match. Complete
GenBank features are retained; internal boundaries are not invented.

A physical SSU link requires a same-specimen alignment at least 97%
identical over at least 80% of the accepted SSU. Marker sequences on
separate contigs can co-occur in a specimen without demonstrably belonging
to the same parasite individual. No new parasite species calls or changes
to the original grades follow from this analysis.

## Outputs and verification

The selected expanded run is
`data/results/waki_comparison/expanded_reference_panel/`. The earlier run
directly under `data/results/waki_comparison/` retains the original
219-reference analysis for provenance. Each directory contains an immutable
script/configuration snapshot and a completion record with input and output
SHA-256 checksums, runtime versions, actual commands, and SLURM job ID.

- `reference_assembly_hits.tsv`: every annotated HSP.
- `long_homology_candidates.tsv`: the length-flagged subset.
- `ssu_assembly_links.tsv`: all accepted-SSU alignments, including failures
  to meet the same-specimen physical-link criteria.
- `report/specimen_marker_summary.tsv`: the 24-row manuscript-facing report.
- `report/all_haplotype_shared_interval_comparisons.tsv`: every eligible
  haplotype scored on the matched contig interval.
- `report/top_marker_fragments.fasta`, `top_marker_pairwise_alignments.fasta`,
  and `top_marker_contigs.fasta`: retained sequence material for review.

The initial `same_interval_comparisons.tsv` table groups GenBank organism
labels and chooses each organism's highest-bit-score HSP first. It is a
retained exploratory diagnostic, superseded for biological comparison by
the report's paper-label grouping and all-haplotype interval calculation.

`tests/test_waki_comparison.py` checks forward and reverse coordinates,
insertions and deletions, invalid coordinates, identical compared
intervals, and alias handling. Seven tests pass in the analysis environment.
Ruff checks and the batch script's shell syntax check also pass. The nine
initial focal alignments were independently checked against the original
GenBank records and Phase-3 contigs, including strand and exact identity
counts. An independent traversal of the expanded panel's raw alignments
also reproduced the Hawaiian/Japanese COI comparisons against the added
*D. barbatus* reference and confirmed no significant Florida match to it.
All hashes in the selected analysis and report completion records verify;
all 24 summary rows have a unique specimen/marker key and recorded host.
