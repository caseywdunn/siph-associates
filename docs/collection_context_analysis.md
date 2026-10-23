# Collection depth, tissue and host context

These descriptive summaries ask which associates occur across sampled depths
and hosts, how the Metamycoplasmataceae bacterial lineage DT-68 varies among
Nanomia collection contexts, and
whether Physalia regional patterns remain visible within shared tissue and host
categories. They retain all **213 libraries**, including zero detections, and use
current v3 grades and reviewed metadata. No additional significance tests or
presence thresholds are introduced by these reporting scripts.

## Inputs and products

[analyze_collection_context.py](../scripts/analyze_collection_context.py) joins
[publication metadata](../figures/cohort/sample_metadata.tsv), accepted microbial,
viral and eukaryotic evidence, physical-flowcell metadata and the focal bacterial
lineage assignments. It checks exact library membership and one-to-one joins.
The source metadata are exported independently of the frozen sequence inputs.
Current source paths and downstream execution are documented in
[expanded reporting](expanded_reporting.md).

Products in [figures/collection_context/](../figures/collection_context/) are:

- `library_context.tsv`: every library, collection evidence, tissue, host,
  region, date, flowcell eligibility and detection indicators.
- `associate_detections.tsv`: accepted evidence with its grades and roles.
- `host_detection_summary.tsv`, `depth_summary.tsv`, `tissue_summary.tsv`:
  denominators and positive specimens for each descriptive grouping.
- `dt68_specimens.tsv`, `dt68_depth_summary.tsv`: the 22 Nanomia and two Resomia
  libraries, including negative specimens.
- `tissue_region_incidence.tsv`, `tissue_region_host_incidence.tsv`: Physalia
  versus recorded tentacles, and the corresponding P. utriculus comparison.
  Single-flowcell subsets are additional descriptions.
- `summary.json`, `provenance.json`: counts, source hashes, command, runtime and
  a copy of the producing script.

[plot_collection_context.py](../figures/plot_collection_context.py) produces
three PDF/PNG figure pairs and their plotting provenance from those tables.

## Collection evidence and interpretation

All 213 libraries now have numeric depth values. The source basis remains
explicit: 150 assigned Physalia surface depths, 40 source-record depths,
14 interval midpoints, one collection-event mean, and eight reported point depths
from Manko2026. One Physalia has a source-record depth of zero. The assigned
surface convention does not become a specimen-specific measurement.

NA19 is 296 m, WS5 is 246 m, and WS6 is the user-specified event mean,
328.589662 m, retaining its 295.8355–350.06003 m bounds. These are resolved values;
older reports that treated them as missing have been superseded. The full
recorded interval determines a depth stratum, even when the point value is a
midpoint or event mean. Among the 62 non-Physalia libraries, 16 are wholly within
0–20 m, 13 have intervals crossing 20 m, and 33 are deeper than 20 m. Depths in
the latter group range from 217 to 3,605.26001 m.

Tissue stays missing for the two unresolved Church OA specimens. The new
Bargmannia libraries have recorded tissue; their retained NCBI collection-source
conflicts are documented in [the study guide](bargmannia.md#collection-metadata-source-conflicts).
The [metadata guide](sample_metadata.md) and cell-level provenance distinguish
reported measurements, reviewed corrections and unresolved identifiers.

Host, locality, tissue, sequencing history and depth overlap. DT-68 comparisons
among Nanomia collection groups do not isolate a depth effect. Positive-specimen
fractions do not standardize detection sensitivity, and tissue comparisons do
not measure bacterial concentration or anatomical location. Genomes, viral
targets and eukaryotic reporting groups also have different counting units.
Use the current generated tables for detection counts rather than earlier
numerical summaries.

## Execution and validation

The [expanded reporting workflow](expanded_reporting.md) connects the tables,
plots and broader overview figures to accepted current inputs. For a bounded
export using already-complete upstream tables, use the
[publication workflow guide](workflow/publication_figures/README.md).
The analysis uses the Python standard library; plotting uses the existing
Snakemake Python environment with Matplotlib and NumPy. Each output records
its actual runtime environment.

Tests in [test_collection_context.py](../tests/test_collection_context.py)
cover identifier joins, point/interval/event-mean depth treatment, accepted
nonredundant evidence and descriptive denominators. Input and output hashes
identify each completed export. Manuscript copying is a later, separate step.
