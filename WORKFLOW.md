# Workflow operation

The paper screens **213 sequencing libraries from four studies** and uses
**209 retained specimens** for biological summaries. Start with
the [paper-to-analysis index](README.md#finding-analyses-from-the-paper) to find
results, or the [workflow index](docs/workflow/README.md) for methods, inputs,
launch commands and generated dependency graphs.

## Select the analysis

The [current cohort guide](docs/expanded_reporting.md) describes read processing,
the completed catalog-v3 mapping and downstream evidence. The current
[`workflow/manuscript_revision.smk`](workflow/manuscript_revision.smk) applies the
exclusion decisions and includes the [matched tissue, host and geographic analyses](docs/coauthor_metadata_analysis.md).
[reporting_sources.json](config/reporting_sources.json) selects the evidence
used by publication exports.

Use each analysis guide's explicit Snakefile, target and configuration. The root
`Snakefile` retains the original 205-library execution; its `all` target checks
readiness and does not reproduce the complete current study. Assembly was
attempted for all 213 libraries. The current
[assembly-completion workflow](docs/assembly_completion.md) has no minimum
retained-pair cutoff.

All commands run from the analysis repository root after the
[environment setup](docs/workflow/README.md#environment-and-inputs). Raw reads,
large databases and completed intermediate files are outside Git; the analysis
guides identify required inputs and the stages that produce them. A bounded
reporting target consumes completed inputs and fails if they are missing.

## Preview and launch

Follow the separate preview and launch commands in the selected guide. In
particular, the [current reporting commands](docs/expanded_reporting.md#preview-and-launch)
use `workflow/manuscript_revision.smk` and the `revision_publication` target.
A dry run describes scheduled work; it does not
establish that the scientific outputs are complete or current.

For documentation-only checks, run:

```bash
./scripts/update_rulegraphs.py --check --require-tracked
```

This checks all configured rule graphs without executing analysis jobs. See the
[graph environment](docs/workflow/graph-environment.md) for the pinned tools and
required data access. Regenerate a changed graph with the same command without
`--check --require-tracked`.

## Monitor and resume cluster work

SLURM is the cluster job scheduler. The stage guides identify the appropriate
launcher and resource allocation; the
[YCRC documentation](https://docs.ycrc.yale.edu/) describes cluster access.
Do not launch overlapping workflows that write the same outputs.

For jobs submitted through `scripts/submit_workflow.sh`, the latest controller
record is `workflow/state/latest_controller.tsv`. Inspect its target,
configuration and commit alongside the scheduler and log:

```bash
column -t workflow/state/latest_controller.tsv
squeue -j CONTROLLER_ID
sacct -j CONTROLLER_ID --format=JobID,JobName,State,Elapsed,ExitCode
tail -n 100 logs/slurm/controller.CONTROLLER_ID.out
```

Other launchers record executions under their analysis result directories, as
specified in their guides. A controller's exit does not necessarily terminate
its independent worker jobs. Before resuming an interrupted run:

1. Let the old controller finish and check that its workers have left both the
   running and pending states. A new controller cannot adopt existing workers.
2. Check the target's validation record and outputs; workers may have completed
   them after the controller stopped.
3. Dry-run the same Snakefile, target and configuration using the original
   launcher environment, then resume that scope only if work remains.
4. If a stale lock remains, confirm that no controller or worker for those
   outputs is active before unlocking that same workflow.

Preserve partial outputs and let the configured `--rerun-incomplete` policy
handle interrupted jobs. Existing launchers use `--rerun-triggers mtime`; a
reviewed change to code or scientific parameters therefore requires explicitly
selecting the affected rules with `-R RULE`. Documentation changes alone do not
require a scientific rerun.

Detailed execution evidence remains in [RUNLOG.md](RUNLOG.md), the analysis
guides' linked run records and per-run provenance. The superseded early gate
instructions are retained in
[the operation archive](dev_docs/workflow-operation-2026-10-23-before-submission.md)
for reconstructing their original execution.
