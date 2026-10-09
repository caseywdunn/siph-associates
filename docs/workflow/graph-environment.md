# Rule-graph rendering environment

The checked-in graphs were generated with the existing McCleary installation.
The documentation tools are separate from the environments used by scientific
rules; changing a graph or README does not change an analysis dependency.

| Component | Pinned generation environment |
| --- | --- |
| Platform | Linux x86_64, RHEL 8.10 |
| Python | 3.11.0 in `/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake` |
| Snakemake | 7.24.0 in that environment |
| PuLP | 2.7.0 in that environment |
| PyYAML | 6.0.3 in that environment |
| Graphviz | `/usr/bin/dot`, RPM `graphviz-2.40.1-45.el8.x86_64` |
| Default `sans` font | Nimbus Sans, RPM `urw-base35-nimbus-sans-fonts-20170801-10.el8.noarch` |
| Fontconfig | 2.13.1-4.el8 |
| Hash ordering | `PYTHONHASHSEED=0` |

Use this installation for both local regeneration and the CI job. Snakemake and
Graphviz versions are enforced by the generator. The font/platform pins matter
for exact SVG geometry; using another renderer or font is not a reason to silently
ignore a freshness failure. The existing `envs/workflow.yaml` describes workflow
setup but differs from the installed Python/PyYAML patch versions above; these
observed versions document the graph build without altering that analysis file.

## CI installation and required data

The GitHub Actions job in
[`rulegraphs.yml`](../../.github/workflows/rulegraphs.yml) runs on an explicitly
labelled **self-hosted Linux runner** (`siph-associates-rulegraphs`) with this
installation and read access to the configured McCleary source/result paths.
The data are too large and partly nonpublic to replace with an invented input set
on a generic GitHub-hosted runner. The workflow uses the checked-out rule/code
files while reading the existing scientific inputs at their configured paths.

The job runs on pushes, pull requests and manual dispatch, uses read-only GitHub
permissions, and executes only the graph freshness command. A fork pull request
must first be reviewed and brought onto a trusted branch; the small hosted gate
fails those requests before any self-hosted job can execute untrusted rule code.

Registering the labelled runner and requiring the `rulegraphs` check in branch
protection are repository-administrator setup steps. This documentation change
does not register a runner, change remote settings or claim a GitHub CI run passed.
The full same-environment check can be run locally using the commands in the
[workflow guide](README.md#regenerate-and-check-the-graphs). If the runner is absent,
the queued CI job is not a successful freshness check.
