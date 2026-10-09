# Workflow operation

For current stage-by-stage methods, input/output descriptions, complete target
coverage and generated rule graphs, use the
[workflow audit guide](docs/workflow/README.md). The operational and historical
gate instructions below remain useful for execution and controller recovery;
their early target list is not the complete set of analyses.

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
- `phase3_benchmark`: run the frozen six-library/12-strategy assembly benchmark,
  including marker, virus, MAG, host-carryover, and resource comparisons. Its
  decision table is the gate for implementing cohort assembly; it is not a
  cohort assembly target.
- `phase3_inputs`: apply the benchmark-selected host route to every library,
  aggregate prepared-pair yields, and validate the objective 175,000-pair
  assembly-eligibility rule before cohort assembly is submitted.
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

The command submits a lightweight seven-day controller on the `week` partition.
The controller launches workers on their rule-specific partitions via
`profiles/slurm/`, so both scheduler and workers survive loss of the interactive
tmux or Codex session. The latest controller ID, target, config, and code commit
are atomically recorded in `workflow/state/latest_controller.tsv`.

On reconnect, inspect the controller, workers, and validated target before doing
anything else:

```bash
column -t workflow/state/latest_controller.tsv
squeue -j CONTROLLER_ID
sacct -j CONTROLLER_ID --format=JobID,JobName,State,Elapsed,ExitCode
tail -n 100 logs/slurm/controller.CONTROLLER_ID.out
test -s tests/work/stages/phase1_smoke.done
```

The profile always uses `--rerun-incomplete`, `--keep-going`, a 180-second
latency wait, no blanket retry, and `--rerun-triggers mtime`. Each output path
belongs to one sample/rule job. Jobs write temporary files before atomic
replacement and emit validated sentinels rather than treating an empty
biological result as success.

### Controller timeout and handoff

Do not cancel a healthy controller merely to move it from `day` to `week`, and
do not start a replacement while workers from the old controller are running or
pending. Generic-cluster workers are independent SLURM jobs; a replacement
Snakemake process cannot adopt them. Because their final outputs do not yet
exist, an overlapping controller can submit duplicate work and cause two jobs
to contend for the same temporary or final paths.

If a controller reaches its walltime, use this procedure:

1. Let the old controller exit. Its already-submitted workers may continue.
2. List workflow workers and wait for every one to leave both the running and
   pending states:

   ```bash
   squeue --me -h -o '%i %j %T' | awk '$2 ~ /^sa\./ && $2 != "sa.controller"'
   ```

3. Check whether the requested target completed while workers drained. For the
   Phase-2 cohort:

   ```bash
   test -s data/results/stages/screen_cohort.done && cat data/results/stages/screen_cohort.done
   ```

4. If the sentinel is absent, dry-run the identical target after the worker
   list is empty:

   ```bash
   /gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/snakemake \
     -n screen_cohort --profile profiles/slurm --configfile config/config.yaml
   ```

5. Submit the identical target with the current committed controller script:

   ```bash
   bash scripts/submit_workflow.sh screen_cohort config/config.yaml
   ```

If that submission reports a stale Snakemake lock, first reconfirm that no old
controller or worker is active, then unlock once and resubmit:

```bash
/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/snakemake \
  --unlock --configfile config/config.yaml
bash scripts/submit_workflow.sh screen_cohort config/config.yaml
```

Do not delete partial outputs during handoff. Atomic publication and
`--rerun-incomplete` make the restart recover accepted work safely.

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

## Phase 3 benchmark gate

The exact sample/strategy matrix is frozen in `config/phase3_benchmark.tsv` and
its parameters, databases, tool prefixes, and resources are isolated in
`config/phase3_benchmark.json`; this avoids invalidating accepted Phase-2
outputs. Validate and launch it with:

```bash
bash scripts/check_phase3.sh --scope static
bash scripts/submit_workflow.sh phase3_benchmark config/config.yaml
```

On completion, require both validators before interpreting the decision table:

```bash
cat data/results/stages/phase3_benchmark.done
python3 scripts/validate_phase3.py --scope benchmark
column -t -s $'\t' data/results/phase3_benchmark/decision_table.tsv | less -S
```

## Phase 3 cohort-input gate

The accepted route choices and eligibility floor are isolated in
`config/phase3_cohort.json`, so they do not invalidate Phase-2 outputs. Validate
the implementation, then prepare and summarize all cohort inputs:

```bash
python3 scripts/validate_phase3.py --scope benchmark
snakemake -n phase3_inputs --profile profiles/slurm --configfile config/config.yaml
bash scripts/submit_workflow.sh phase3_inputs config/config.yaml
```

After completion, require `data/results/stages/phase3_inputs.done` to report
PASS and inspect `data/results/phase3_cohort/input_eligibility.tsv`. Freeze that
table before submitting any cohort assemblies.

## Phase 3 cohort assembly

`workflow/rules/phase3_assembly.smk` assembles each eligible library and
recovers its rRNA markers, viruses, and MAG bins. The tool commands are
visible in the rule `shell:` blocks. Membership is frozen from the accepted
eligibility table, and settings and resources live in
`config/phase3_assembly.json`; neither file invalidates accepted outputs.

| Step | Rule | Main output under `data/results/phase3_cohort/` |
|---|---|---|
| Assembly (≥1 kb contigs, IDs prefixed with the library) | `assemble_contigs_megahit` | `assemblies/<sample>.fasta` |
| rRNA markers (bac/arc/euk) | `find_rrna_barrnap` | `markers/<sample>.gff` |
| Viral contigs | `identify_viruses_genomad` | `viruses/<sample>.fna`, `.tsv` |
| Viral quality | `assess_viruses_checkv` | `checkv/<sample>.tsv` |
| Coverage binning (BAM kept only in scratch) | `bin_contigs_metabat2` | `bins/<sample>.tar.gz`, `.depth.tsv` |
| Bin quality | `assess_bins_checkm2` | `checkm2/<sample>.tsv`, `.status` |
| Cohort table, including below-floor exclusions | `summarize_assemblies` | `assembly_summary.tsv` |
| Membership and accounting check | `validate_assemblies` | `../stages/phase3_assembly.done` |

Freeze membership, check, dry-run, and launch from the repository root:

```bash
python3 scripts/freeze_phase3_assembly.py          # production; --fixture for tests
bash scripts/check_phase3.sh
snakemake phase3_assembly -n --profile profiles/slurm --configfile config/config.yaml
bash scripts/submit_workflow.sh phase3_assembly config/config.yaml
```

After completion, `python3 scripts/validate_phase3.py --scope assembly` checks
every expected output and assembly checksum independently. This stage reports
QC values without applying thresholds: MAG/vOTU quality rules, pooled GTDB-Tk,
and dereplication are locked in the following stage.

## Phase 3 catalogs

`workflow/rules/phase3_catalog.smk` builds the cohort MAG and vOTU catalogs
from the accepted assembly products. It applies the locked rules in
`config/phase3_catalog.json`, with evidence and reasoning in
`docs/phase3_qc_decisions.md`.

| Step | Rule | Output under `data/results/phase3_catalog/` |
|---|---|---|
| MAGs ≥50% complete, <10% contaminated | `select_mags` | `mags/candidates.tsv` |
| GTDB r220 placement (GTDB-Tk 2.4.1) | `prepare_gtdbtk_reference`, `classify_mags_gtdbtk` | `mags/gtdbtk.*.summary.tsv` |
| Pairwise ANI | `compare_mags_skani` | `mags/candidates.skani.tsv` |
| Prokaryotic MAGs, species at 95% ANI / 50% AF | `build_mag_catalog` | `mags/mag_catalog.tsv`, `mags/species_representatives/` |
| Viral inclusion (hallmark and ≥5 kb or ≥50% complete) | `select_viruses` | `viruses/included.tsv` |
| Host alignment and endogenous-candidate flags | `align_viruses_host_minimap2`, `flag_endogenous_viruses` | `viruses/virus_classification.tsv` |
| vOTUs (95% ANI, 85% AF) | `cluster_votus_checkv`, `build_votu_catalog` | `votus/{associate,endogenous_candidate}.votus.tsv` |
| Rule and accounting check | `validate_catalog` | `../stages/phase3_catalog_<catalog_version>.done` |

```bash
snakemake phase3_catalog -n --profile profiles/slurm --configfile config/config.yaml
bash scripts/submit_workflow.sh phase3_catalog config/config.yaml
```

## Phase 4 catalogs

`workflow/rules/phase4_catalog.smk` freezes the version-v1 mapping catalogs
under `config/phase4_catalog.json`, with reasoning in
`docs/phase4_catalog_decisions.md`. The accepted Phase-2 nomination table is
checked against its SHA-256 instead of being declared as an input, so the
known trim-provenance mtime cascade cannot reopen the screens.

| Output under `data/results/phase4_catalog/v1/` | Contents |
|---|---|
| `bacterial_catalog.fna` (+ bwa index, `.fai`) | Catalog genomes and decoys; contigs named `<catalog_id>\|<contig>` |
| `bacterial_catalog.manifest.tsv` | Every nominated, MAG, and decoy genome, with status, substitute, taxonomy, CheckM2, evidence, and checksums |
| `viral_associate.fna`, `viral_endogenous_candidate.fna` (+ manifests) | Frozen vOTU representatives with sequence checksums |
| `crispr/host_links.tsv` | vOTU–host links from full-length spacer matches (≤1 mismatch) |

```bash
snakemake phase4_catalog -n --profile profiles/slurm --configfile config/config.yaml
bash scripts/submit_workflow.sh phase4_catalog config/config.yaml
```

## Environments and provenance

Pinned portable definitions are under `envs/`, and the Yale profile activates
the managed Conda environments. The legacy phyloFlash child process uses its
tested cluster prefix because its EMIRGE dependency requires Python 2, while
the workflow wrapper remains in pinned Python 3.
Every substantive job writes a per-rule log, benchmark, software/parameter JSON,
and checksum-linked provenance. Each run snapshots its config, sample table, and
Phase-0 manifest freeze under its persistent result root.
