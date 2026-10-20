# Review candidate genomes before full-cohort mapping

The previous Thermoproteota assignment was a Perkinsela-related eukaryotic bin.
GTDB placement and CheckM2 estimates alone did not establish prokaryotic identity;
the [domain audit](mag_domain_audit.md) documents the evidence. Every prospective
MAG catalog member, including previously retained bins, needs further review.

[Preparation](assembly_completion.md#prepare-the-combined-cohort-for-qc) preserves
sequence identities and source assemblies. Known eukaryotic bins stay excluded
from the prokaryotic catalog and remain available as positive controls.

## Required evidence

1. **Identity and provenance:** link each candidate FASTA to its library and bin
   archive with checksums; verify unique contig IDs and distinguish replacement
   assemblies. Reuse old taxonomy only for identical sequences.
2. **Gene-domain evidence:** extend DIAMOND/Swiss-Prot searches to every candidate
   and use nr to adjudicate discordant evidence. Report bacterial, archaeal, eukaryotic, viral, ambiguous
   and unmatched genes. Retain database versions, controls and near-best-hit
   ties. Check translated DNA windows where prokaryotic gene prediction may
   miss eukaryotic genes. No hit does not establish prokaryotic identity.
3. **Contig consistency:** compare gene-domain evidence by contig with GC,
   coding density, coverage and length. Investigate mixed-domain bins,
   host-derived sequences and classifications supported by only a few contigs.
   GC differences alone do not justify removal.
4. **Marker support:** examine GTDB marker counts, alignment support, warnings
   and agreement with broader gene evidence. Conserved marker hits alone do
   not establish domain. Review every proposed archaeal assignment and conflict.
5. **Genome quality:** apply existing completeness and contamination criteria
   after establishing prokaryotic identity. Record retained, excluded and
   unresolved decisions with evidence, then dereplicate accepted genomes.
6. **Viruses and references:** review viral hallmarks, CheckV quality,
   host/endogenous-sequence evidence (including the new Halistemma reference),
   and nominated reference-genome identities and provenance.

Preparation labels gene-domain and contig-taxonomy QC as pending. Assembly
completion and candidate extraction do not complete these reviews. The final
review must produce per-candidate decisions and an immutable inventory of
approved mapping targets with sequence checksums and inclusion evidence.

No preparation workflow includes full-cohort mapping rules. The explicit hold
in `config/post_expansion_remapping.json` cannot be cleared by preparation.

## October 20 review checkpoint

The completed broad screen covers all 131 candidates, using both predicted
proteins and translated DNA windows against the frozen Swiss-Prot database.
Fifteen domain controls passed. The
[review table](../figures/catalog_qc/domain_review.tsv) combines domain evidence,
marker counts, coding density and assembly coverage. It remains an evidence
table, not an approved catalog.

There are 104 bacterial-majority and 27 eukaryotic-majority candidates. The latter
include the four known Perkinsela-related controls, 20 previously unclassified
bins, and all three new Bargmannia candidates. New GTDB classification supports
the additional Physalia bin from YPM-IZ-111016 as Tyloplasma; the three new
Bargmannia bins have no bacterial or archaeal markers.

The nr follow-up contains 6,406 proteins across all 131 candidates: every
eukaryotic/archaeal Swiss-Prot assignment in a bacterial-majority candidate, all
genes on 15 contigs with at least three eukaryotic assignments exceeding their
bacterial assignments, and the 20 longest eukaryotic-hit proteins from each
eukaryotic-majority candidate. These counts select evidence for review; they are
not automatic exclusion thresholds. This supplementary nr search uses DIAMOND
`--sensitive`; the whole-candidate Swiss-Prot searches use `--very-sensitive`.
The original four-bin nr audit remains available separately.

Completed preparation also includes all three host-reference comparisons for
750 initially selected viral sequences, retained geNomad gene annotations and
CheckV boundary evidence, candidate genome ANI comparisons, and quality/taxonomy
assessment of 34 newly nominated references (22 GTDB and 12 OceanDNA). The eight
Bargmannia screen outputs enter the same nomination parser as the original 205
libraries; original nominations and reference QC are reused with provenance.
The negative-control reference pool remains fixed.

Before approving mapping targets, finish nr adjudication, resolve discordant
contigs, review the viral cellular flanks and hallmark evidence, combine old and
new reference QC, and curate/dereplicate the accepted sequences. Reference merge
requires compact bacterial and archaeal GTDB summary tables; the new full GTDB
output is retained under `phase4_catalog/references/new/gtdb/`. Final catalog
decisions, target checksums, and index preparation are still pending. No
full-cohort mapping has been started.

## Run and inspect the bounded workflows

From the repository root, use the Snakemake 7.24.0 environment described in the
[workflow guide](workflow/README.md#environment-and-inputs). The launchers below
capture source/configuration snapshots, logs and exit status. They are independent
entry points consuming the staged products; run the broad domain screen before
its nr follow-up. None imports cohort read-mapping rules.

| Purpose | Preview without executing | Cluster launcher |
| --- | --- | --- |
| All-candidate domain evidence | `snakemake domain_screen -s workflow/catalog_qc.smk --cores 1 --nolock --rerun-triggers mtime --dry-run` | `sbatch scripts/catalog_qc.sbatch` |
| nr adjudication | `snakemake nr_review -s workflow/catalog_nr_review.smk --cores 1 --nolock --rerun-triggers mtime --dry-run` | `sbatch scripts/catalog_nr_review.sbatch` |
| Viral catalog evidence and new-bin taxonomy | `snakemake catalog_evidence -s workflow/post_expansion_catalog.smk --cores 1 --nolock --rerun-triggers mtime --dry-run` | `sbatch scripts/post_expansion_catalog.sbatch` |
| New nominated reference QC | `snakemake new_reference_evidence -s workflow/expanded_references.smk --cores 1 --nolock --rerun-triggers mtime --dry-run` | `sbatch scripts/expanded_references.sbatch` |
| Viral genes and boundaries | `snakemake viral_evidence -s workflow/viral_catalog_qc.smk --cores 1 --nolock --rerun-triggers mtime --dry-run` | `sbatch scripts/viral_catalog_qc.sbatch` |

Broad-domain and nr outputs are in `data/results/catalog_qc/2026-10-20/`, with
execution records for completed broad-domain jobs 13633299/13633697 and nr
controller 13633804. Other outputs are in `data/results/post_expansion/2026-10-20/`,
with successful controller records 13633504 (catalog evidence), 13633990
(references) and 13634130 (viral genes/boundaries). Large sequence/search outputs
and execution archives are retained on the cluster outside Git.

The graphs below are generated from those exact targets and tracked configurations.
Their freshness is checked by the existing rule-graph CI job. Regenerate with
`./scripts/update_rulegraphs.py catalog_qc catalog_nr_review post_expansion_catalog expanded_references viral_catalog_qc`.

![All-candidate domain evidence](workflow/graphs/catalog_qc.svg)

![nr follow-up](workflow/graphs/catalog_nr_review.svg)

![Initial catalog evidence](workflow/graphs/post_expansion_catalog.svg)

![New nominated reference evidence](workflow/graphs/expanded_references.svg)

![Viral gene and boundary evidence](workflow/graphs/viral_catalog_qc.svg)
