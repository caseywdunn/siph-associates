# Analysis workflow audit guide

Start here to follow the analyses from the original sequencing libraries to
the reported figures. Each of the 19 rule files has a descriptive README below,
with inputs, products, scientific steps, exact rule names, preview/launch commands
and a generated dependency graph. The guides describe the current implementation;
the linked decision records explain the accepted scientific choices.

The source workflow, analysis configurations and accepted results are unchanged
by this documentation. Viewing or regenerating the graphs does not run analyses.
Production `all` means **input readiness only**, not the full paper workflow.

## Follow the data

| Stage / rule file | What to audit | Guide |
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
See the [project manifest guide](../../README.md#the-manifest) and
[`data/metadata/`](../../data/metadata) for input inventory and provenance.
Fresh execution requires access to those files and the configured databases;
the small smoke test also draws reads from real source libraries. Some later
stages deliberately require previously accepted, checksum-identified products.
No graph command downloads or invents missing inputs.

Each stage guide gives a copyable dry run and a separate cluster launch command.
Most targets can construct missing upstream products. The publication guide
instead uses an explicit list of permitted reporting rules, matching the
existing controller's supported bounded-execution mechanism.

## Read the graphs

The 16 distinct scope graphs are generated directly with `snakemake --rulegraph`
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

The guides cover all current `workflow/rules/*.smk` files. They do not imply that
every reported output is built by Snakemake:

- Manifest construction, metadata enrichment, surface-depth scoring and the
  assembly-membership freeze use separate scripts and reviewed source tables.
- The direct Waki marker comparison has its own
  [methods and reproduction guide](../waki_sequence_comparison.md).
- The current study overview and broad host-distribution figures use the explicit
  [overview plotting command](../../figures/README.md#biological-overview-figures).
- Selected assets are copied to `../manuscript_siph_associates`; manuscript
  compilation and editorial changes happen in that separate repository.

The [run log](../../RUNLOG.md), [figure provenance](../../figures/README.md) and
stage decision records identify completed, accepted outputs. Older v1 outputs
and historical plan sections are retained; current mapping and inference use v2.

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
