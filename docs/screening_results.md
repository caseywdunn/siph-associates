# Reporting the initial screening results

The manuscript's screening comparison summarizes existing, accepted tables.
It distinguishes classifier nominations, genus corroboration, catalog inclusion,
competitive-mapping presence and host-aware SSU assignments. No classifier,
assembly, alignment, statistical model or presence grading is rerun.

Generate the reporting bundle with Python 3.10 or later, using only the standard
library:

```bash
python3 scripts/summarize_screening_results.py
```

Outputs in `figures/screening/` comprise generated LaTeX values, a JSON summary,
method-label counts, retained sylph nomination/library pairs with their existing
catalog-v2 grades, SSU control comparisons, source database and corroboration
summaries, input/output hashes, and a snapshot of the producing script.
The bundle also contains accepted prokaryotic detections partitioned by
same-library sylph nomination, nomination in another library, or MAG-only
recovery; the accepted eukaryotic evidence with assembly sources; and broad
Mycoplasma-related screen support for libraries with accepted focal MAGs.

The source tables are:

- `data/results/aggregation/cohort/{library_qc,candidate_nominations}.tsv`:
  accepted library scope and all returned Bracken, sylph and phyloFlash labels.
  Bracken is the abundance re-estimation stage following Kraken2; it is not an
  independent classifier replicate. Taxonomic labels are not accepted richness.
- `data/results/phase4_catalog/nomination_evidence/`: source-database and
  same-genus corroboration summaries. Only the GTDB subset has this particular
  cross-taxonomy corroboration assessment.
- `data/results/phase4_catalog/v2/bacterial_catalog.manifest.tsv` and
  `data/results/phase6_analysis/v2/grades/bacterial_grades.tsv`: final retained
  representatives and the existing same-library nomination and presence grades.
  Excluded references and multiple nominated genomes merged into one catalog
  species are not separate tests in the mapping comparison.
- `data/results/phase6_analysis/v2/sensitivity/leave_one_study_out.tsv`:
  independent retained accounting check for genomes contributed only by MAG
  recovery, without a sylph nomination anywhere in the cohort.
- `data/results/eukaryote_gate/read_evidence.tsv`: paired-read assignments from
  the initial phyloFlash mappings, after correcting the historical taxonomy
  extraction issue; this is the pre-competitive-remapping comparison.
- `data/results/eukaryote_gate/read_lineages.tsv`,
  `calibration/negative_controls.tsv`, `eukaryote_grades.tsv`, and
  `config/eukaryote_gate.json`: subsequent host-aware assignments, background
  calibration, corrected nonredundant reporting and fixed thresholds.
- `config/phase6_analysis.json`: the existing minimum aligned assembly bases
  used to summarize contig support for accepted sylph nominations.

The script checks unique library/catalog coverage, agreement of Bracken row
counts with accepted QC, agreement of raw sylph counts with the retained source
summaries, MAG-only accounting and the remapped controls' calibration values.
It checks input hashes again before recording successful completion.

## Interpretation

The strongest evidence of false-positive screening assignments concerns the
widespread initial insect SSU signal that diminishes below the presence
threshold after host-aware competitive remapping. This supports misassignment
of conserved host reads to non-host references as an important source of
spurious apparent associates. The original decision record is
[`eukaryote_gate_decisions.md`](eukaryote_gate_decisions.md).

Genome nominations that fail the subsequent breadth/read criteria remain
unsupported at those criteria, rather than established false positives. Low
abundance, sequence divergence, competing reference genomes and differences in
host handling can reduce later support. Genus corroboration also depends on
database coverage and naming. The study does not provide a separate measured
false-positive rate for Kraken2, Bracken or sylph, nor a truth-labelled benchmark
of all nominated organisms. De novo MAGs without sylph nominations illustrate
the complementary sensitivity of assembly.

The overlap tables distinguish two kinds of additional recovery: a genome
absent from all sylph nominations, and a genome nominated elsewhere but not
in the particular library where mapping detects it. Broad Mycoplasma-related
support uses genus names containing `mycoplasm` in Bracken (at least ten
estimated reads) and SSU lineages containing that string in phyloFlash (at
least one read). These are the same permissive support levels used for genus
corroboration, not new presence thresholds or matches to specific MAGs.

Eukaryotic assembly overlap counts only nonredundant accepted detections.
An assembly-only confirmation can still have weak read evidence: the
YPM-IZ-111009 kinetoplastid example has supporting reads below the accepted
threshold, with its confirming SSU from the general assembly rather than
phyloFlash. The Discussion connects the different perspectives supplied by
read-based taxonomic profiles, genome-reference nominations, ribosomal
markers, assembly and competitive validation.

The manuscript's bounded transfer is:

```bash
cd ../manuscript_siph_associates
python3 scripts/sync_assets.py --screening-only --dry-run
python3 scripts/sync_assets.py --screening-only
```

This copies only the new reporting bundle and writes its own checksum manifest,
preserving other figure and metadata exports. LaTeX reads the copied values;
building the manuscript does not access the analysis repository or run analyses.
