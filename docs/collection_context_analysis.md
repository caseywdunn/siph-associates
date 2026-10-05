# Collection depth, tissue and host context

## Scope and scientific questions

The 2026-10-04 metadata enrichment permits three exploratory descriptions:
which organisms were recovered across sampled depths and hosts; whether DT-68
varies among Nanomia collection contexts; and whether the previously observed
Physalia regional pattern is visible within a shared tissue and host category.
These analyses retain all 205 libraries. They use accepted primary detections,
introduce no new presence thresholds or significance tests, and do not rerun
sequence analyses. Source base revision is
`7cac98ca4c82f1e8c76c8bd714f6510b753416b1`; new code and outputs are identified
by exact hashes and source snapshots, rather than the base revision alone.

## Inputs, outputs and interpretation

`scripts/analyze_collection_context.py` joins `manifest.csv` to accepted
exports in `figures/cohort/`: bacterial, viral and eukaryote detections,
physical-flowcell metadata and the verified Mycoplasmatales species mapping.
It validates a one-to-one library join and keeps libraries with zero accepted
detections. All outputs are in `figures/collection_context/`:

- `library_context.tsv`: all 205 libraries, depth provenance and intervals,
  tissue, host, region, collection date, flowcell eligibility and detections.
- `associate_detections.tsv`: accepted evidence with original grades and roles.
- `host_detection_summary.tsv`, `depth_summary.tsv`, `tissue_summary.tsv`:
  specimen denominators and detections for each descriptive grouping.
- `dt68_specimens.tsv`, `dt68_depth_summary.tsv`: all 22 Nanomia and two Resomia.
- `tissue_region_incidence.tsv`: all Physalia versus recorded tentacles;
  `tissue_region_host_incidence.tsv`: the same comparison within P. utriculus.
  Single-flowcell subsets are additional descriptions, not new tests.
- `summary.json`, `provenance.json`: selected counts, methods, input/output
  hashes, actual command, runtime and producing source snapshot.

`figures/plot_collection_context.py` consumes these retained tables and produces
three vector PDFs, PNG previews and `plot_manifest.json`. Genomes, viral targets
and eukaryotic groups have different counting units. Positive-specimen
fractions do not standardize detection sensitivity. Missing tissue stays unknown;
depth intervals have no invented midpoint. All 150 assigned Physalia surface
depths remain distinct from the one recorded Physalia depth of zero. NA19,
WS5 and WS6 remain unresolved. Collection years for NA34–NA36 retain the
existing values with their source disagreement documented in the metadata audit.

## Selected observations

Among non-Physalia, 29 specimens were collected at 0–20 m, 22 at 217–1349 m,
and three have unresolved depths. DT-68 occurs in all four N. septata with
accepted deeper depths, none of the three shallow specimens, and all three
with unresolved depths. The positive specimens came from California and
the negative specimens from Washington. All ten have the same primary study
and recorded siphosome tissue, but locality, collection history and sequencing
batch differ. None of the four accepted deeper specimens is single-flowcell
eligible. They meet mapped-presence criteria but not the high-confidence
assembly criterion; WS7 retains subthreshold assembly evidence.

The southwestern-versus-central Pacific contrast for MAGSP0005/0007/0010
persists within Physalia tentacles (18/30, 10/30, 18/30 versus 0/15, 1/15,
1/15) and P. utriculus tentacles (3/8, 4/8, 4/8 versus 0/15, 1/15, 1/15).
The eastern Indian Ocean has no recorded tentacle samples. Western Indian
Ocean tentacle counts remain high, but all five specimens share one flowcell.
These results refine hypotheses about associations; they do not identify
independent effects of depth, tissue or geography.

## Reproduction and validation

See [the figure build guide](../figures/README.md) for the bounded Snakemake
target, environment and exact preview/run commands. The scripts are also
individually executable from the repository root. Analysis requires Python
3.10+ standard library; the recorded run used Python 3.11.4. Plotting used
Matplotlib 3.7.1 and NumPy 1.24.2. Each producer records its actual runtime.

Nine scientific tests in `tests/test_collection_context.py` cover consequential
identifier joins, depth treatment, evidence selection and subset denominators.
They passed, together with Ruff and a fresh temporary replay whose ten
scientific outputs were byte-identical. Independent review recomputed detection
counts and metadata joins from the original accepted inputs and verified every
regional host-subset denominator and numerator. Existing primary sequence and
statistical outputs remain unchanged. Figure checks and the manuscript copy
manifest complete the handoff; the manuscript compiles from copied assets.
