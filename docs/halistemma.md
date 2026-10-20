# Halistemma rubrum host depletion

**Current status (2026-10-20):** both depletion mappings completed. The original
driver stopped at its assembly-count gate. That cutoff has been withdrawn;
both libraries now proceed through [assembly completion](assembly_completion.md)
using all retained pairs. Use that guide for the current assembly and QC route.
The workflow and launch commands below document the original depletion run.

Two existing Ahuja2024 libraries are reprocessed using the newly public
[GCA_987347065.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_987347065.1/)
reference, jhHalRubr1.1. Their stable sample identifiers are unchanged:

| Library | Voucher | Locality | Previously quality-filtered pairs |
| --- | --- | --- | ---: |
| Ahuja2024__CWD12 | YPM-IZ-35038 | Villefranche, Mediterranean | 100,871,052 |
| Ahuja2024__CWD17 | YPM-IZ-35040 | Northwestern Atlantic | 165,797,762 |

Both libraries are below the 200-million-raw-pair limit. Their original
25-million-pair full-read assemblies remain available for comparison. The
new run asks whether removing host reads allows better associate assembly
from the complete libraries. It does not repeat the completed fish screen.

## Workflow

The [workflow](../workflow/halistemma.smk) reuses the established trimming,
BWA depletion, MEGAHIT, barrnap, geNomad, CheckV, MetaBAT2 and CheckM2 rules.
Results are isolated under `data/results/halistemma/2026-10-19/`.

1. Download the pinned assembly and check NCBI's compressed-file MD5, accession,
   species, sequence IDs and lengths. Include the nuclear scaffolds and
   mitochondrial sequence listed in the assembly report. Build the BWA index.
2. Regenerate full trimmed reads with the same fastp settings. Check original
   raw-file hashes and read counts against the completed fish-screen records.
3. Draw exactly one million paired reads per library with BBTools seed
   20261019 and map them using the same BWA implementation as full depletion.
   Report mapped reads and the fraction of pairs retained with both mates
   unmapped. Validate paired IDs and BAM/FASTQ accounting.
4. Require at least half the pilot pairs to have a mapped mate in each
   library before automatically proceeding. This is an operational check that
   host removal is useful, not a biological detection threshold. A lower value
   stops the workflow for review; it does not exclude the specimen.
5. Map all trimmed pairs and retain pairs with both mates unmapped, using the
   existing `samtools fastq -f 12` convention. Assemble all retained pairs when
   the established 175,000-pair eligibility minimum is met, then recover
   markers, viral sequences and bins and assess their quality.

Full trimmed and host-depleted reads are retained for the combined analysis.
The launcher archives source, metadata, parameters and execution logs. Reference
size is 4.38 Gb; indexing and mapping request 64 GB RAM. Mapping benchmarks and
pilot summaries provide the evidence for subsequent resource adjustments.

From the repository root, using the established Snakemake environment:

```bash
XDG_CACHE_HOME="$PWD/.cache" snakemake --snakefile workflow/halistemma.smk \
  halistemma --cores 16 --dry-run
mkdir -p logs/halistemma
sbatch scripts/halistemma/launch.sbatch
```

The bounded target `halistemma_pilot` stops after the two pilot mapping reports
and their check. The default `halistemma` target continues through assembly QC
only if those checks pass. Neither target includes pooled catalog construction
or cohort-wide competitive remapping.

![Halistemma workflow](workflow/graphs/halistemma.svg)

## Hold on combined remapping

**Do not remap the full cohort until additional MAG/viral QC, catalog review,
and explicit mapping authorization are complete.** This extends the user's
2026-10-19 hold with their 2026-10-20 instruction, recorded
in [post_expansion_remapping.json](../config/post_expansion_remapping.json).

Review the eight new Bargmannia libraries and both replacement Halistemma
assembly inputs together, including domain assignments and possible host
carryover. Select one assembly input per Halistemma library for the final
catalog; the original assemblies remain as comparison baselines, not duplicate
inputs with overlapping contig names. Then build and freeze the pooled microbial and viral catalogs,
retaining the MAGSP0026 domain exclusion. The [QC plan](catalog_qc.md) precedes
authorization for one full-cohort remapping of **213** libraries. Halistemma
changes processing for two existing samples; it does not add samples.

The current 205-library manuscript remains the reporting baseline until that
combined analysis passes review. Neither of the active per-library workflows
automatically starts remapping.
