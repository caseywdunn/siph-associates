# Crustacean COI pilot

This analysis tests whether COI searches recover crustacean sequences in
host-depleted *Nanomia septata* libraries. It compares exact nucleotide seeds
with translated protein searches on the same read pairs. It is separate from the
accepted manuscript analyses and does not change their inputs or identifiers.

## Inputs and references

[The pilot manifest](../config/crustacean_coi_pilot.tsv) includes all ten
*N. septata* libraries and every retained host-depleted pair: 2,152,939 pairs in
total, with no additional subsampling. These are the existing
`data/results/phase3_cohort/inputs/{sample}_R[12].fastq.gz` files. They retain pairs
with both mates unmapped to the host after quality filtering and the original
200-million-pair input cap. Neither host mapping nor earlier analyses are rerun.
The [broader screening plan](prey_marker_screen.md) defines the staged expansion.

The [download specification](../config/crustacean_coi_sources.json) includes
[MetaZooGene](https://metazoogene.org/atlas/) COI barcodes for copepods, krill,
decapods, amphipods, mysids, lophogastrids and ostracods. Annotated COI sequences
from the [RefSeq mitochondrial release](https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/)
extend coverage across the gene. Other RefSeq COI sequences and MetaZooGene
cnidarian COI barcodes provide competitors during nucleotide identification.
Downloaded URLs, dates, HTTP headers and SHA-256 hashes are recorded in
`data/results/crustacean_coi/reference/2026-10-17/downloads.json`.

Source headers are retained in the reference metadata. The downloaded MetaZooGene
FASTA headers identify `MZGdb_v2025-m10-01`; the website's newer database label is
not treated as the version of these files. Barcodes can occupy different reading
frames and orientations: DIAMOND against annotated crustacean COI proteins identifies their
frames, followed by translation with invertebrate mitochondrial code 5. Accepted
barcodes require at least 50 aligned amino acids, 50% protein identity, 80% query
coverage and no internal stop codon. Only unambiguous sequences of 180–2,400 nt
enter this pilot panel. Exclusions and group coverage are tabulated.

## Search and interpretation

1. BBDuk recruits pairs with an exact 25-base match on either strand to the
   DUST-masked target COI panel. Both mates are retained.
2. DIAMOND searches each mate in translation using code 5, sensitive mode,
   minimum ORF length 20, at least 60% query coverage and E-value at most 0.001.
   The protein panel contains the same target records as the nucleotide panel.
3. Both mates of the union of candidates are aligned with BLASTN against the
   target and competitor COI sequences. The extraction checks pair synchronization
   and the expected input count.
4. Provisional classification requires at least 90 aligned nucleotides and 80%
   identity. References within 1% of the best qualifying bit score are retained.
   Target/non-target ties and discordant mates are ambiguous. A result list
   truncated while scores are still near-best is also ambiguous.

These are pilot comparison criteria, not a validated species-identification or
prey-presence threshold. A COI match does not by itself establish predation.
Performance depends on reference coverage and read depth;
an absent pilot match does not establish absence from the specimen.

Twelve paired synthetic controls exercise exact and reverse-complement recovery,
synonymous substitutions, siphonophore host sequences, fish/insect COI, random
sequence and low complexity. Real-read processing depends on these checks passing.
Per-library `summary.json`, `classified_pairs.tsv` and benchmark files report
method overlap, competitive assignments and runtime separately.

## Environment and execution

The existing McCleary environments supply Python 3.11, Biopython 1.81, Snakemake
7.24.0 and BBTools 39.01. The launcher loads BLAST+ 2.15.0; DIAMOND 2.1.24 and
SeqKit 2.8.0 (installed under module label 2.8.1) are selected by absolute paths.
Paths and dated output roots are in [the configuration](../config/crustacean_coi.json).
For another installation, rebuild these dependencies with conda-forge/bioconda
and update those paths before running.

From the repository root:

```bash
export PATH=/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:$PATH
export XDG_CACHE_HOME="$PWD/.cache"

# Preview the isolated workflow only.
snakemake --snakefile workflow/crustacean_coi.smk pilot \
  --cores 8 --dry-run --rerun-triggers mtime

# Run references, controls and all ten complete retained libraries on YCGA.
mkdir -p logs/crustacean_coi
sbatch scripts/crustacean_coi/batch.sh pilot

# Optional smaller scopes use the same rules.
sbatch scripts/crustacean_coi/batch.sh references
sbatch scripts/crustacean_coi/batch.sh controls
```

Run only one scope at a time against these output directories. Existing successful
products are reused. Logs and execution-specific source snapshots, input file
statistics, versions, exit status, commands and report copies are retained under
`data/results/crustacean_coi/executions/{SLURM_JOB_ID}/`.
Bulk reference and read files remain outside Git.

![COI pilot dependencies](workflow/graphs/crustacean_coi.svg)

This graph is generated from `workflow/crustacean_coi.smk`, target `pilot`; repeated
sample/mate jobs collapse into one node per rule. It does not invoke the main
analysis workflow. Controls and real libraries use the same search rules, so
collapsing them creates an apparent loop at the control gate; the actual jobs
form an acyclic dependency graph. Generate/check it without running any analysis:

```bash
./scripts/update_rulegraphs.py crustacean_coi
./scripts/update_rulegraphs.py crustacean_coi --check
```

The existing rule-graph CI check covers this entry. Detailed execution history is
in [the pilot run notes](../dev_docs/crustacean-coi/pilot.md).
