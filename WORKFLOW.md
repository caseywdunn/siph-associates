# Workflow operation

The root `Snakefile` is the publication workflow. It reads only the normalized
sample table and locked configuration; no rule reads the exploratory repository.
Production products are written under `data/results/`, while regenerable trimmed
reads are written under `/vast/palmer/scratch/dunn/cwd7/siph_associates/`.

## Targets

- `phase1_ready`: validate the production sample table and snapshot configuration.
- `phase1_smoke`: run the three-study, three-host-route fixture through trimming,
  Kraken2/Bracken, sylph, phyloFlash, host handling, and output validation.
- `screen_pilot`: run universal trimming and all three screens for the eight
  full-library cases frozen in `config/phase2_pilot.tsv`, then validate and
  aggregate them without Phase-3 host handling.
- `screen_cohort`: run the same universal screens, validation, library-QC
  aggregation, and candidate-nomination aggregation for all 205 libraries. It
  must not be submitted until the full-library pilot resource gate passes.
- `all`: resolves to `phase1_ready` in production and `phase1_smoke` in the fixture
  configuration.

The fixture is deterministically drawn from 50,000 paired reads from these
predeclared real libraries:

| Study | Fixture | Host route |
|---|---|---|
| Church 2025 | `Church2025__FM-16644` | `P_physalis` reference |
| Ahuja 2024 | `Ahuja2024__CWD1` | reference-free/Kraken classified-out |
| Ahuja 2026 | `Ahuja2026__NA19` | `N_septata` reference |

Fixture reference mapping uses deterministic 1 Mb excerpts solely to exercise
indexing and both-unmapped-pair handling. Production configuration points to the
complete audited references.

## Validation before submission

```bash
bash scripts/check_phase1.sh
```

This rebuilds both sample tables, rechecks the Phase-0 freeze, runs Snakemake
lint, dry-runs the 205-library and fixture DAGs, and writes
`data/metadata/phase1_validation.txt`. After the smoke finishes, run:

```bash
bash scripts/check_phase1.sh --require-smoke
python3 scripts/freeze_phase1.py
```

## Persistent SLURM execution

Submit the fixture controller from the repository root:

```bash
bash scripts/submit_workflow.sh phase1_smoke tests/config.fixture.yaml
```

The command submits a lightweight 24-hour controller. The controller launches
workers via `profiles/slurm/`, so both scheduler and workers survive loss of the
interactive tmux allocation. The latest controller ID, target, config, and code
commit are atomically recorded in `workflow/state/latest_controller.tsv`.

On reconnect, inspect the controller, workers, and validated target before doing
anything else:

```bash
column -t workflow/state/latest_controller.tsv
squeue -j CONTROLLER_ID
sacct -j CONTROLLER_ID --format=JobID,JobName,State,Elapsed,ExitCode
tail -n 100 logs/slurm/controller.CONTROLLER_ID.out
test -s tests/work/stages/phase1_smoke.done
```

If the controller times out, submit the identical target again. The profile
always uses `--rerun-incomplete`, `--keep-going`, a 180-second latency wait, no
blanket retry, and `--rerun-triggers mtime`. Each output path belongs to one
sample/rule job. Jobs write temporary files before atomic replacement and emit
validated sentinels rather than treating an empty biological result as success.

Trimmed reads are normal scratch outputs and downstream rules consume them with
`ancient()`. In production Phase 2 they are marked temporary and are removed
only after Kraken2/Bracken, sylph, phyloFlash, and paired-count validation have
all succeeded; this bounds residency below the shared scratch capacity. A
purged trimmed file is regenerated only if a pending downstream target needs
it, and regenerating scratch reads does not invalidate already-complete
persistent products. Capped libraries are limited to one concurrent trim
because exact seeded sampling temporarily stages their concatenated raw mates.
When code or locked parameters change, record and force the affected rules
explicitly with `-R RULE`, because production uses mtime-only rerun triggers.

## Phase 2 gate

Run the static and fixture checks before the full-library pilot:

```bash
bash scripts/check_phase2.sh --scope static
bash scripts/submit_workflow.sh screen_pilot config/config.yaml
```

The pilot covers eight declared libraries and reuses accepted persistent
outputs in the cohort. After it passes, use rule benchmarks and `sacct` maxima
to lock resources with margin, rerun `scripts/check_phase2.sh --scope pilot`,
and commit the resource decision before submitting `screen_cohort`.

## Environments and provenance

Pinned portable definitions are under `envs/`. The Yale profile uses the tested
cluster installations declared in `config/config.yaml`; rule code first honors
an activated pinned environment and otherwise resolves the declared prefix.
Every substantive job writes a per-rule log, benchmark, software/parameter JSON,
and checksum-linked provenance. Each run snapshots its config, sample table, and
Phase-0 manifest freeze under its persistent result root.
