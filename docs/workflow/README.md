# Analysis methods and rule graphs

For the current 213-library execution, start with
[current cohort and downstream reporting](../expanded_reporting.md).
The original stage guides below explain the shared algorithms and historical
205-library gates. The accepted expansion uses those same rules through these
bounded entry points:

| Current operation | Guide | Workflow |
| --- | --- | --- |
| Fourth-study reads, screens and assemblies | [Bargmannia](../bargmannia.md) | `workflow/bargmannia.smk` |
| Halistemma host depletion | [Halistemma](../halistemma.md) | `workflow/halistemma.smk` |
| Assemble every retained input | [Assembly completion](../assembly_completion.md) | `workflow/assembly_completion.smk` |
| Genome/viral QC and catalog freeze | [Catalog QC](../catalog_qc.md) | QC, curation and freeze entry points linked in that guide |
| Completed all-library mapping | [Remapping](../cohort_remapping.md) | `workflow/cohort_remapping.smk` |
| Grades, statistics, SSU evidence and exports | [Reporting](../expanded_reporting.md) | `workflow/expanded_reporting.smk` |


Follow DNA from the original sequencing libraries to the paper's figures and
supporting tables. The [paper-to-analysis index](../../README.md#finding-analyses-from-the-paper)
connects Results headings to these methods. The guides below cover all 21 rule
files: each explains the biological purpose, inputs, outputs and interpretation,
and places a generated rule graph beside its preview and run commands.

Read screening and assembly are complementary branches. All libraries can be
compared with shared microbial sequences, including those without an assembly.
A rule graph shows how computational steps depend on one another; it does not
show that a run completed or imply that every sample yielded a genome.

## Follow the data

| Stage / rule file | What it explains | Methods and graph |
| --- | --- | --- |
| `common.smk` | Sample membership, read routes, configuration and default target | [Shared decisions](common/README.md) |
| `preprocess.smk` | Read cap, adapter/quality filtering and small test inputs | [Read preparation](preprocess/README.md) |
| `screens.smk` | Kraken/Bracken, sylph and phyloFlash nominations | [Universal screens](screens/README.md) |
| `host.smk` | Host indexing and the original smoke-test handling route | [Initial host handling](host/README.md) |
| `validation.smk` | Configuration snapshot, readiness and smoke checks | [Initial validation](validation/README.md) |
| `phase2.smk` | Pilot/cohort screen aggregation and quality accounting | [Screen aggregation](phase2/README.md) |
| `phase3.smk` | Six-library comparison of assembly-input strategies | [Discovery benchmark](phase3/README.md) |
| `phase3_cohort.smk` | Accepted input strategy and assembly eligibility | [Cohort inputs](phase3_cohort/README.md) |
| `phase3_assembly.smk` | Contigs, ribosomal markers, viral sequences and genome bins | [Assembly and recovery](phase3_assembly/README.md) |
| `phase3_catalog.smk` | Genome quality, species clustering and viral clustering | [Recovered catalogs](phase3_catalog/README.md) |
| `phase4_catalog.smk` | Public references, decoys, shared catalogs and CRISPR links | [Mapping catalogs](phase4_catalog/README.md) |
| `phase5_mapping.smk` | All-library competitive mapping and continuous coverage | [Cohort mapping](phase5_mapping/README.md) |
| `phase6_analysis.smk` | Presence grades, geographic/batch inference and controls | [Grades and statistics](phase6_analysis/README.md) |
| `eukaryote_gate.smk` | Conservative non-host SSU evidence and reporting | [Eukaryotic associates](eukaryote_gate/README.md) |
| `mycoplasmatales_phylogeny.smk` | Relationships inferred from conserved bacterial proteins | [Genome phylogeny](mycoplasmatales_phylogeny/README.md) |
| `mycoplasmatales_gene_content.smk` | Predicted functions and limits of absence claims | [Gene content](mycoplasmatales_gene_content/README.md) |
| `cnidarian_16s.smk` | Similar markers in other sampled cnidarians | [Cnidarian context](cnidarian_16s/README.md) |
| `cassiopea_16s.smk` | Independent jellyfish tissue/environment comparison | [Cassiopea context](cassiopea_16s/README.md) |
| `publication_figures.smk` | Bounded exports from accepted results | [Figures and source tables](publication_figures/README.md) |
| `fish_12s_search.smk` | Four fish-marker searches and competitive identification | [Fish 12S workflow](../prey_marker_screen.md) |
| `host_depletion_report.smk` | Per-library input, removal and retained-read counts | [Host-depletion statistics](../host_depletion_report.md) |

A sequencing library is the unit processed by the workflow. A MAG is a
metagenome-assembled genome recovered by grouping assembled sequences; a vOTU is
a cluster of related viral sequences; SSU denotes small-subunit ribosomal RNA.
Candidate nominations, competitive read detections and assembly-supported
detections are different evidence levels. Genome discovery can miss microbes
that are absent from the reference catalog even though every library is mapped
against every catalog member. These distinctions matter when comparing host groups.

## Environment and inputs

All commands in the guides run from the **analysis repository root**, not from
the README's directory. On the existing McCleary installation:

```bash
cd /gpfs/ycga/work/dunn/cwd7/siph_associates
export PATH=/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin:$PATH
export XDG_CACHE_HOME="$PWD/.cache"
unset SNAKEMAKE_PROFILE SIPH_ASSOCIATES_ALLOWED_RULES
snakemake --version
dot -V
```

The workflow uses Snakemake 7.24.0. Individual rules refer to environments in
[`envs/`](../../envs) and to installed tool prefixes in their configuration files.
The existing cluster launcher uses [`profiles/slurm/`](../../profiles/slurm) and
the resource overrides in
[`workflow_controller.sbatch`](../../scripts/workflow_controller.sbatch).
Those settings control scheduling; they add no scientific branches to these
graphs. SLURM is the cluster job scheduler; consult the
[YCRC cluster documentation](https://docs.ycrc.yale.edu/) for access and setup.
The operational details and recovery procedure remain in
[`WORKFLOW.md`](../../WORKFLOW.md).

Original FASTQs and large databases are external to Git. Their locations are
recorded in [`config/samples.tsv`](../../config/samples.tsv),
[`config/config.yaml`](../../config/config.yaml) and the stage JSON files.
See the [sample and input overview](../../README.md#data-and-interpretation) and
[`data/metadata/`](../../data/metadata) for input inventory and provenance.
Fresh execution requires access to those files and the configured databases;
the small smoke test also draws reads from real source libraries. Some later
stages deliberately require previously accepted, checksum-identified products.
No graph command downloads or invents missing inputs.

Each stage guide gives a copyable dry run and a separate cluster launch command.
Most targets can construct missing upstream products. The publication guide
instead uses an explicit list of permitted reporting rules, matching the
existing controller's supported bounded-execution mechanism.

### Rebuild the original sample manifest (historical reconstruction)

The retained 205-library manifest is a frozen input to the original analyses.
Do not rebuild or replace it for publication-metadata edits or expanded reporting.
For a deliberate reconstruction of that original input version,
use [build_manifest.py](../../scripts/build_manifest.py), which combines the
source-study tables, reads sequencing headers from the original FASTQs, and
applies the [reviewed metadata overlay](../../data/metadata/sample_metadata_updates.tsv).
This requires access to the recorded sequencing paths. Resource inventory also
requires the configured reference databases. From the repository root in the
environment above:

```bash
python3 scripts/build_manifest.py
python3 scripts/inventory_phase0_resources.py
python3 scripts/summarize_phase0.py
python3 scripts/validate_manifest.py
```

These commands write `manifest.csv` and inventory, provenance, summary and
validation tables in `data/metadata/`. They do not run sequence analyses and
have no dry-run mode. Check the regenerated inputs against the recorded freeze:

```bash
python3 scripts/validate_manifest.py --freeze data/metadata/manifest.freeze.sha256
```

If intentionally establishing a new accepted input version after reviewing the
differences, `python3 scripts/freeze_manifest.py` records new checksums. Replacing
the freeze does not demonstrate reproduction of the original inputs. The
[metadata source review](../../data/metadata/enrichment/2026-10-04/README.md)
documents the retained collection records and unresolved discrepancies.

## Read the graphs

The scope graphs are generated directly with `snakemake --rulegraph`
and rendered with Graphviz. The preprocessing, host and initial-validation guides
share the same three-library fixture graph; screening and Phase 2 share the
pilot/cohort graph. `common.smk` contains no rules, so its guide shows the readiness
target selected by its configuration. Repeated library jobs collapse into one
node per rule. Arrows point from a producer to a consumer; colors distinguish
rules and have no biological meaning. Open an SVG directly to zoom into larger
graphs. The corresponding `.dot` files preserve Snakemake's exact graph output.

Graphs describe declared dependencies, not execution success or all files read
inside scripts. They do not display input datasets as nodes. In particular:

- Frozen screen nominations enter reference nomination through a checksum-checked
  parameter, and accepted phyloFlash archives enter eukaryote extraction through
  provenance. These handoffs do not create upstream edges.
- The benchmark decision and assembly-membership freeze are reviewed transfers
  into configuration, not automatically inferred by downstream rules.
- Several selection, extraction and plotting scripts read additional files through
  directory parameters or fixed paths. Each affected guide identifies that limit.
- Genome selection is shared by phylogeny and gene-content analyses; the final
  inferred tree is not an input to gene prediction.

The dependency graphs retain these properties of the existing workflow. An absent
edge does not mean the underlying scientific evidence was independent. Consult
the source rules, scripts and accepted provenance when checking an individual claim.

## Analyses outside these rule files

The bounded [expanded-catalog QC workflows](../catalog_qc.md) have separate
entry points and generated graphs. They review all prospective genomes, new
reference nominations and viral boundaries before the now-completed full-cohort remapping.

The guides cover all current `workflow/rules/*.smk` files. They do not imply that
every reported output is built by Snakemake:

- Manifest construction, metadata enrichment, surface-depth scoring and the
  assembly-membership freeze use separate scripts and reviewed source tables.
- The direct Waki marker comparison has its own
  [methods and reproduction guide](../waki_sequence_comparison.md).
- The broad host-distribution and eukaryote figures are now produced by
  `expanded_overview` in the [reporting entry point](../expanded_reporting.md).
- Selected assets are copied to the separate manuscript repository. Compilation
  and the explanatory methods drawing (Figure 1) are maintained there; they
  are not outputs of these rule files.

The [figure guide](../../figures/README.md) identifies selected outputs and
their source tables. Current mapping and inference use catalog v3. Detailed
execution evidence is retained in the [run log](../../RUNLOG.md) and stage records.

## Regenerate and check the graphs

The shared executable is outside the scientific dependency graph. It reads
[`graphs.json`](graphs.json), which records only targets, configuration and the
publication rule restriction; graph edges always come from Snakemake.

```bash
# With the environment above active, generate every documentation graph.
./scripts/update_rulegraphs.py

# Regenerate into temporary files and compare; do not replace existing graphs.
./scripts/update_rulegraphs.py --check

# Select one scope, or preview its normal execution plan without running jobs.
./scripts/update_rulegraphs.py mapping --check
./scripts/update_rulegraphs.py mapping --dry-run

# CI additionally rejects graphs missing from Git's tracked files.
./scripts/update_rulegraphs.py --check --require-tracked
```

The generator fixes `PYTHONHASHSEED=0`, checks the Snakemake/Graphviz versions,
keeps diagnostics separate, and stages all requested outputs before replacement.
It neither executes rule bodies nor invokes the cluster launcher. Normal dry runs
can report jobs that would be rerun under current file timestamps; those are plans,
not permission or a requirement to execute them for this documentation change.

The renderer/platform pins and CI prerequisites are in
[`graph-environment.md`](graph-environment.md). CI uses the same generator and
fails on stale or untracked DOT/SVG files. Include regenerated artifacts with
future workflow changes. No automatic commit or branch-protection change is made.
The generation and review record is in
[`dev_docs/workflow-documentation.md`](../../dev_docs/workflow-documentation.md).
