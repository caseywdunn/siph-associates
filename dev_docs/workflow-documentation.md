# Workflow documentation verification, 2026-10-09

The user requested descriptive READMEs and generated rule graphs for every rule
file, following the Dunn Lab
[rule-graph guidance](https://github.com/caseywdunn/dunnlab_code/blob/main/skills/dunnlab-workflow-design/references/rulegraphs.md),
without changes requiring analysis reruns. The linked reference was fetched and
read on 2026-10-09. The scientific source base is
`24336ab0343d436079d5abe3ef8d6e15a0fdc10e`.

## Added documentation

[`docs/workflow/README.md`](../docs/workflow/README.md) indexes 19 per-file guides.
Together they describe all 133 rules in `workflow/rules/*.smk`; `common.smk` has
no rules and documents shared configuration instead. Sixteen distinct target
scopes have generated DOT/SVG pairs. Shared fixture/screen scopes are reused by
the corresponding guides rather than drawn independently. Each guide uses the
same targets, configuration and (for publication) rule restriction for its
launch instructions, dry run and graph.

The source Snakefile, rule files, configurations, environments, analysis scripts
and scientific outputs were not edited. README/operation-guide navigation was
updated, including the root README's stale Phase-3-active status. The graph
generator and CI configuration are outside the scientific dependency graph.
No cluster job, analysis rule, metadata freeze or manuscript build was executed.

## Verification

- `scripts/update_rulegraphs.py` generated all 16 scopes using Snakemake 7.24.0
  and Graphviz 2.40.1. A second complete pass with `--check` reproduced every
  DOT and SVG byte-for-byte. No DOT normalization or hand-maintained edges are used.
  A narrow `.gitattributes` exemption preserves Snakemake's closing-line spaces
  in generated DOT files while allowing Git whitespace checks elsewhere.
- `scripts/update_rulegraphs.py --dry-run` completed for every documented scope.
  The bounded publication scope was up to date. Some unbounded scopes plan work
  because the current snapshot inputs are newer or temporary reads are missing;
  none of those planned jobs was executed. Documentation does not resolve or
  invalidate historical execution state.
- Every rule name in each source file is described in its matching guide and
  appears in that guide's SVG. The common-helper exception is labelled explicitly.
- Generated SVG XML and local Markdown link targets were checked. The large
  Phase-6 graph, the fixture graph and the bounded publication graph were also
  inspected visually. Larger SVGs should be opened directly for readable labels.
- The new Python generator passed Ruff formatting and E/F/I checks. Ruff was
  installed under `/tmp` only; no analysis environment was changed.
- A deliberately altered SVG was rejected by `--check` without being overwritten;
  restoring the original made the check pass. The CI tracked-file guard also
  correctly rejected the new, uncommitted graph files. CI YAML parsed with the
  intended push/pull-request/manual triggers and the freshness command.

Local temporary build/preview/check logs are `/tmp/siph-rulegraphs-build.log`,
`/tmp/siph-rulegraphs-preview.log` and `/tmp/siph-rulegraphs-check.log`. They are
diagnostic aids, not scientific provenance or required inputs. The stable
regeneration specification is `docs/workflow/graphs.json`; observed renderer and
font versions are in `docs/workflow/graph-environment.md`.

## Interpretation boundaries exposed by the audit

The guides explicitly describe handoffs that are not declared Snakemake edges:
accepted screen nominations, phyloFlash archive provenance, directory-based
sequence extraction/selection, and the reviewed assembly-membership freeze.
They also distinguish the historical reference-free Kraken route from the
accepted cohort assembly subset and the fuller read input used for mapping.

The current overview plots and Waki marker comparison are outside the rule files
and link to their existing reproduction records. The Cassiopea guide records
that Illumina and PacBio branches have different alignment-coverage criteria;
this is a description of existing code, not a scientific-method change.

## CI setup boundary

The freshness workflow is configured for pushes, pull requests and manual runs
using a labelled self-hosted runner with the recorded renderer and access to
the existing inputs. Generic hosted CI cannot reconstruct these graphs from
the repository alone. The job rejects missing/untracked or stale graph files.
The initial documentation pass left the new files uncommitted for review. The
user subsequently requested a commit and push; the tracked-file condition requires
the graph files to be included in that commit. No GitHub run is claimed here.

An administrator must provide the `siph-associates-rulegraphs` runner and make
the `rulegraphs` check required if merge enforcement is desired. No remote
runner or branch-protection setting was created in this task.
