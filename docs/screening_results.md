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
catalog-v3 grades, SSU control comparisons, source database and corroboration
summaries, input/output hashes, and a snapshot of the producing script.
The bundle also contains accepted prokaryotic detections partitioned by
same-library sylph nomination, nomination in another library, or MAG-only
recovery; the accepted eukaryotic evidence with assembly sources; and broad
Mycoplasma-related screen support for libraries with accepted focal MAGs.

The current paths are selected explicitly by
[reporting_sources.json](../config/reporting_sources.json):

- Expanded nominations under `data/results/post_expansion/2026-10-20/aggregation/cohort/`
  combine the original 205 and eight Bargmannia libraries using the shared parser.
- `data/results/post_expansion/2026-10-21/reporting/screens/` contains the
  213-library screening inventory and refreshed source-database/corroboration
  summaries. Bracken row counts are checked against the original screen inventory
  and the new study's completed Bracken reports.
- [The reviewed v3 catalog](../figures/catalog_qc/catalog_v3/) and current
  `phase6_analysis/v3/grades/bacterial_grades.tsv` provide catalog membership,
  same-library nominations and accepted detections. Excluded references and
  genomes merged into a representative are not separate mapping tests.
- Current `phase6_analysis/v3/sensitivity/leave_one_study_out.tsv` checks
  nomination dependence and MAG-only accounting across all four studies.
- Current `eukaryote_gate/read_evidence.tsv` combines the original initial
  phyloFlash assignments with those from the eight new archives. Competitive
  `read_lineages.tsv`, negative-control summaries and final grades use the
  expanded host-aware reference.
- [eukaryote_gate.json](../config/eukaryote_gate.json) and
  [phase6_analysis.json](../config/phase6_analysis.json) retain the locked
  read-level and assembly-support thresholds.

The `expanded_screening_summary` rule connects these products in the
[expanded reporting workflow](expanded_reporting.md). It can be previewed there;
the standalone script above only reads completed tables.

The script checks unique library/catalog coverage, agreement of Bracken row
counts with accepted QC, agreement of raw sylph counts with the retained source
summaries, MAG-only accounting and the remapped controls' calibration values.
It checks input hashes again before recording successful completion.

## Interpretation

In the original 205-library analysis, the small initial insect SSU signal
diminished below the presence threshold after host-aware competitive remapping.
This supports misassignment of conserved host reads to non-host references as
an important source of spurious apparent associates. The expanded cohort adds
one exception to the original background result: formalin-preserved Bargmannia
elongata T676-SS3 has 288 insect and 71 plant pairs at ≥97% identity, above the
locked 37-pair threshold. Its primary grades remain unchanged, with a specimen
flag and a descriptive sensitivity summary omitting it in the
[expanded comparison](../figures/expanded_review/README.md). The original decision record is
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

The expanded cohort has one archaeal reference detection: `BAC00450`
(`GCA_013390375.1`, Nitrosopelagicus) in Bargmannia stenotes D666-D10. It was
nominated by sylph and passes the primary read-mapping rule (13.1901% breadth,
1,083 retained reads), with no passing same-library assembly support. It remains
below the 20% breadth sensitivity threshold. The screening summary counts this
under archaeal references, not archaeal MAGs. No archaeal MAG is recovered;
the previously misassigned MAGSP0026 remains excluded.

The manuscript's bounded transfer is:

```bash
cd ../manuscript_siph_associates
python3 scripts/sync_assets.py --screening-only --dry-run
python3 scripts/sync_assets.py --screening-only
```

This copies only the new reporting bundle and writes its own checksum manifest,
preserving other figure and metadata exports. LaTeX reads the copied values;
building the manuscript does not access the analysis repository or run analyses.
