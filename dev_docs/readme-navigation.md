# Manuscript navigation revision, 2026-10-11

The README now addresses a biologist reading the paper who wants to understand
its methods and locate the supporting evidence. This revision follows the
[specified Dunn Lab documentation guidance](https://github.com/caseywdunn/dunnlab_code/blob/de3c03d07356390c6b17e3e8aea388e97bf9a70b/skills/dunnlab-workflow-design/references/documentation.md#connect-the-manuscript-to-the-repository)
and the [hypx README example](https://github.com/caseywdunn/hypx/blob/main/README.md).
Both were retrieved through the GitHub API and read on 2026-10-11.
The retrieved text SHA-256 values were:

- Guidance: `7ada9eeabd6cc25baa6d47ca1bb4af10f952ca148b87b1b2a1a0187bca4e3228`.
- Example: `8dad46297eb02996f277b65db80eca7fca4b706a55a2ac12e2b5a3b7a03e022b`.

The analysis base revision was `00a629cd917440a97a8a5562b70e6837a5115d83`. The current
manuscript source was inspected at base revision `0cb4e12710757cf4985b70ae4d4856ff53495401`,
including its uncommitted Dalila revisions. The manuscript title, Results headings,
figure order, supplementary tables and selected asset transfer record were the
basis of the navigation; the earlier analysis report was not used as the paper.

## Scope

- Rewrote the root README around the study design, three chronological source
  studies, the nine manuscript Results headings, and nine major analysis steps.
  Each Results row links to methods and retained evidence. Each analysis step
  links to its guide and generated graph.
- Clarified the parallel read-screening and assembly routes, host-reference use,
  descriptive versus inferential sampling, and limits on interpretation and
  reproduction. Missing public manuscript/archive identifiers remain placeholders.
- Updated the workflow and figure indexes, publication guide, collection-context
  input description and links in the operation and rules guides. Manifest rebuild
  instructions now have their canonical home in the workflow index.
- Removed outdated figure-selection and execution-history prose from the reader
  entry points. Existing historical artifacts, execution records and graph
  verification records were retained. No scientific assets were deleted.

## Verification

- Counted distinct library and specimen identifiers in `manifest.csv`: 205 of each.
  Primary-study assignments are Ahuja 2024: 33, Church 2025: 151, Ahuja 2026: 21,
  in chronological table order. Reused samples are counted once.
- Confirmed 151 Physalia and 22 Nanomia libraries; host-reference routes comprise
  151 Physalia and ten N. septata libraries, leaving 44 without host references.
- Compared all nine Results headings and their order to the current manuscript
  sources, expanding the separate Mycoplasmatales Results input.
- Matched the nine repository-generated figure links to the manuscript's figure
  order and selected source paths in its asset manifest. Figure 1 is explicitly
  identified as a manuscript-owned explanatory drawing.
- Checked 168 local Markdown links and ten heading anchors in the seven edited
  guides, table column consistency, and tracked targets for the root README.
- Confirmed that all 19 rule-file guides exist, embed graphs, and provide separate
  dry-run and launch commands. All 16 graph scopes have tracked nonempty DOT/SVG
  pairs; every SVG parses as XML.
- Confirmed unchanged hashes for 154 workflow, configuration, reporting-source
  and graph files captured before editing. Existing work on figures, metadata
  and reporting scripts was preserved. `git diff --check` passed.

No analysis, cluster job, figure build, manifest rebuild or graph regeneration was
run for this documentation change. The graph and launch commands retain the
[previous generation and dry-run verification](workflow-documentation.md).
This revision checked graph availability and unchanged content; it does not claim
a new Snakemake freshness run or execution of the full workflow. The one-off
link/count checks were run against both the staged and installed documentation.
