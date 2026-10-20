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
2. **Gene-domain evidence:** extend the existing DIAMOND/Swiss-Prot and nr audit
   to every candidate. Report bacterial, archaeal, eukaryotic, viral, ambiguous
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
