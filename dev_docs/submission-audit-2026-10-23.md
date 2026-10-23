# Submission correspondence audit, 2026-10-23

This audit compares the current paper and S1 Appendix with the 213-library
reporting selected by `config/reporting_sources.json`. It checks reporting and
reader access; it does not repeat sequence processing, assembly or mapping.

## Corrections

- The README now follows all 12 Results headings in their current order, with
  main Figs 1–4/Table 1 and appendix Figs A–H/Tables A–D.
- Current guides use the reviewed Nanomia identifications, viral classifications
  and detection counts, collection metadata and completed reporting status.
  The external cnidarian comparison's high-identity abundance maximum is 0.38%.
- Reporting instructions give the prerequisite order for expanded statistics,
  reviewed catalogs/grades, expanded publication products and final reporting.
- Historical tree-label quality estimates and Waki alignment host names are
  distinguished from current genome-quality and specimen metadata. Original
  sequence comparisons and frozen evidence remain unchanged.
- The workflow index covers all 22 shared rule files. The regenerated expanded
  publication graph removes two obsolete upstream edges. An unreferenced old
  collection-context graph was removed; the maintained graph covers that scope.
- CI now makes completed results available at the relative path required by
  downstream rules. Superseded operational and reader guides are preserved as
  explicitly dated records.

## Verification

All 257 selected manuscript assets match their analysis sources, and 238 recorded
source/output hashes pass. All 72 generated display commands reproduce exactly.
Independent counts from the detection and sample tables agree with the paper;
manual checks cover regional statistics, sensitivity comparisons, genome and
marker comparisons, fish and human screening, and Methods thresholds/versions.
The manuscript corrections are the old Nanomia depth-comparison wording and the
version notation for lme4 (2.0-6) and vegan (2.7-5).

All 711 local links and anchors in 63 public reader documents pass, and the 62
public file/directory paths cited by the paper exist. The README title and all
12 Results headings match the paper. All 43 DOT/SVG rule-graph pairs reproduce
in an isolated checkout with the recorded tools and completed-result inputs.
The same checkout fails before those inputs are made available, confirming the
CI repair. No scientific jobs were executed.

The main, illustrated-review and appendix PDFs compile with resolved citations
and cross-references and no LaTeX warnings or box errors; both Word files also
rebuild. Rendered checks include the corrected depth paragraph, Nanomia labels,
main evidence table and appendix figures/tables. Publication source and Word
links do not point to the private manuscript repository.

## Remaining limits

OA-190531-1 and OA-Bluebottle still lack verified public run accessions. A complete
public deposit of reconstructed genomes and markers remains pending; the public
supporting export contains selected trematode sequences. The archived v0.5.0
code release predates current reporting. These gaps remain explicit in the paper
and reader documentation and prevent calling the submission package complete.

Remote GitHub CI execution has not been verified. Graph construction was checked
locally; an attempted all-scope Snakemake dry run stopped because sandbox access
made external outputs appear protected. This is not a completed workflow rerun.
