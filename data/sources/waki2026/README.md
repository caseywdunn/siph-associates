# Waki et al. 2026 reference sequences

References for direct comparisons with the trematodes reported in
[Waki et al., Journal of Helminthology 100:e4](https://doi.org/10.1017/S0022149X25100989).
This archive contains **108 newly deposited sequences: 26 28S, 11 amplicons
containing ITS2, and 71 COI sequences**. It retains fish-derived adults and
larvae and larvae from cnidarians and a chaetognath. These are reference
sequences, not 108 distinct parasites or species.

`references.fasta` uses the unique, versioned GenBank accession as each
record ID. `references.gb` preserves the complete retrieved GenBank records.
`reference_metadata.tsv` links every sequence to the paper's taxon label,
GenBank organism, host, isolate, developmental stage, collection date,
location, marker, length, and sequence checksum. Hosts and stages come from
each accession's source feature; a species having another host in the paper
does not assign that host to every sequence. `reference_features.json`
preserves all feature annotations, their original Biopython location strings,
and zero-based, end-exclusive coordinates.

The selection uses all 108 nucleotide records linked by NCBI ELink to the
article's PubMed record, PMID 41424079, independently verifying completeness
against the printed lists. These comprise 98 unique LC889 accessions in the
13 sheets of the original supplementary workbook, eight accessions explicitly
listed in the paper's DNA-marker sections (LC889172, LC889180, LC889181,
LC889188, and LC889268–LC889271), and two additional records omitted from
both printed lists: LC889189.1 (Lepocreadiidae sp. 1, 28S, a *Physalia*
metacercaria) and LC889237.1 (*Dinurus barbatus*, COI, a fish adult).
All 108 linked GenBank records cite the exact article title.
`pubmed_nuccore_links.json` preserves the publication-to-sequence links,
and metadata distinguishes the two database-only accessions.
Supplementary accession/label rows, including older references
not selected for this study-only FASTA, are retained in
`supplement_accessions.tsv`. Six supplementary rows have an asterisk rather
than a deposited accession; these are retained in
`supplement_unaccessioned.tsv` without inventing additional records.

The workbook is preserved as `S0022149X25100989sup001.xlsx`.
`retrieval_manifest.json` records its publisher URL, the NCBI EFetch query,
retrieval or cached-file timestamps, checksums, execution command, Python and
Biopython versions, and sequence counts. Source downloads are cached; an
explicitly removed download is retrieved again on the next run. The
supplement's expected checksum is pinned so a changed workbook requires
review before its accession assignments are accepted.

From the analysis repository root, reproduce the archive with the existing
Snakemake Python environment (**Python 3.11.0, Biopython 1.81** in this
retrieval) and the system `curl`:

```bash
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python scripts/waki_comparison/retrieve_waki_references.py
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python scripts/waki_comparison/retrieve_context_references.py
```

This script retrieves, parses, and validates the references; it does not run
alignments. Sequences are checked for nonempty nucleotide content, unique
versioned IDs, exact requested/retrieved accession correspondence, article
attribution, and consistency between supplementary and GenBank marker labels.
An independent parse with openpyxl 3.1.5 matched every supplementary
accession row; counts and the workbook checksum are recorded in
`supplement_parser_validation.json`. FASTA, GenBank, and metadata sequence
identities, lengths and hashes were also cross-checked for all 221 records.

The study reference set contains no 18S sequences. **ITS2-labelled records
include 5.8S and 28S flanks, and GenBank does not annotate the internal
boundaries.** A short match to one of these records therefore does not by
itself demonstrate a match within ITS2. Assess the aligned region before
using such a match for species identification. Closely matching a reference
also requires considering coverage, alternatives, and physical linkage to
the specimen's other parasite evidence.

## External species comparisons

The separate `context_references.fasta` and `context_references.gb` contain
all **113 uniquely accessioned external references** used in the supplement:
43 labelled 28S, 22 ITS2, and 48 COI. These include *Dinurus longisinus*
AY222202.1, *D. euthynni* OP458333.1, earlier *D. tornatus* sequences, and
the other *Prodistomum* species available in Waki's comparisons. They are
not specimens newly sequenced for Waki's study. `context_metadata.tsv` and
`context_features.json` retain the same accession/marker/host links and
feature annotations; `context_retrieval_manifest.json` records selection,
versions, checksums, and source query. The main 108-reference files are
unchanged when retrieving this comparison set.

Paper taxon labels and current GenBank organisms are both retained, including
previous unnamed labels; their differences are not silently resolved.
Supplement marker labels can describe only one part of a longer deposited
sequence. Inspect the feature annotations when interpreting a short match.

The paper's *Prodistomum* sp. 1 accession paragraph repeats LC889249–LC889250
for both adults and larvae. The supplement and GenBank instead identify
LC889248.1 as a larva from *Urashimea globosa*, with LC889249.1 and
LC889250.1 from fish-derived adults. This archive retains those
accession-specific GenBank annotations rather than propagating the repeated
range. LC889252.1 (*P. orientale* from *Clytia* sp.) is also included from
the supplement although absent from the main article's accession paragraph.

The original 106-reference selection, derived from the paper and supplement
before the independent PubMed-link check, is preserved byte-for-byte under
`initial_106_references/`, including its retrieval script and provenance.
It was used only in the initial exploratory comparison. The expanded
108-reference archive includes the omitted *D. barbatus* COI record and is
the appropriate input for final comparisons; omission of that reference
must not be interpreted as an absence of a matching marker in our samples.
