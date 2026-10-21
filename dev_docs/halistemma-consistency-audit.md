# Halistemma processing audit, 2026-10-21

Casey Dunn requested equal treatment of Halistemma and the other host-mapped
libraries, followed by commit and push to `main`. This audit follows the
[completed fish reruns](fish-12s/2026-10-21-halistemma-host-depleted.md).

## Findings and corrections

The selected sequence analyses already use the correct depleted inputs. Two
bookkeeping/entry-point inconsistencies remained:

- `scripts/reporting_manifest.py` inherited `host_reference=none` for CWD12 and
  CWD17 from the frozen original manifest. It now joins current host routes from
  `config/post_expansion_samples.tsv` through the reviewed stable analysis IDs.
  Original manifests, read paths, specimen handles and sample selection remain
  unchanged. Mismatched or duplicate join keys fail explicitly.
- The Halistemma launcher still defaulted to its historical assembly target,
  whose 175,000-pair gate could reject CWD12 on reproduction. Its default now
  stops at `halistemma_depletion`; the documented next step is the existing
  cutoff-free assembly-completion workflow. Historical targets/settings remain
  available for reproducing their original executions. No depletion command,
  reference, parameter or selected assembly changed.

The original fish-input manifest producer's help now explicitly identifies its
frozen 161-library Physalia/Nanomia scope. Halistemma uses the same fish search
rules with its separate two-library input manifest. Current documentation and
the generated Halistemma graph reflect these boundaries.

## Evidence examined

| Stage | Verified handling |
| --- | --- |
| Host depletion | Same shared `phase3_prepare_cohort_input.py`, BWA `mem`, samtools 1.9, and `fastq -f 12` as Physalia and N. septata; both mates must be unmapped |
| Assembly and recovery | CWD12/CWD17 MEGAHIT logs use the complete retained FASTQs; shared marker, geNomad, CheckV, MetaBAT2 and CheckM2 rules; no read-count cutoff |
| Selected catalog inputs | Both staged assembly paths resolve to completion outputs, and their SHA-256 values match `assembly_sources.tsv`; original full-read assemblies are excluded from selection |
| Host contamination | `H_rubrum` is present in the selected catalog/curation configurations and completed candidate-contig and viral host-alignment evidence |
| Microbial and viral mapping | Both executed BWA command logs use the depleted paths; sequence counts equal twice the retained pair counts for each library and both catalogs; shared filtering/deduplication/coverage rules |
| Assembly support and SSU | Current reporting consumes the staged depleted assemblies; competitive SSU reads come from the original full-trimmed phyloFlash screen, as for every host |
| Initial screens | Kraken2, Bracken, sylph and phyloFlash use full quality-filtered inputs before depletion for all hosts |
| Fish 12S | Both completed replacements use every retained pair and the shared search/competitive-classification methods |
| Mapping summary | Both are included in the 163-library table; individual-read counts equal twice the corresponding pair counts |

CWD12 contributes 131,703 retained pairs and CWD17 contributes 2,352,421 to
assembly, microbial mapping, viral mapping and fish screening. Their selected
assembly hashes are respectively
`49e265c9164b130a60b619ea00283d9a8353e6cbefac52162a6aac4a639d1db9` and
`9196e3ca3ca4a5c5d1c548dc9285f0e66e70e79f90afb8fe4c2c0ccc2df00a4a`.

The evidence paths are recorded by `config/post_expansion_mapping_inputs.tsv`,
`data/results/post_expansion/2026-10-20/assembly_sources.tsv`, and the logs under
`data/results/assembly_completion/2026-10-20/logs/phase3_cohort/` and
`data/results/post_expansion/2026-10-20/remapping/logs/phase5_mapping/`.
Historical v2 read-handling controls retain their original fixed subsets; they
are not current per-library analyses and are not relabelled as v3 results.

## Validation and execution

The 37 relevant tests passed (fish classification/report selection, current-route
reporting joins, Halistemma inputs, host-depletion reporting, remapping inputs,
and expanded reporting). Ruff formatting/lint and launcher syntax checks pass.

The reporting dry run scheduled 20 reporting/statistics steps and no sequence
analysis. Allocation **13691813** completed all 20 with exit status zero. Source,
logs and execution status are under
`data/results/post_expansion/2026-10-21/reporting/executions/13691813/`.
That record also preserves `pre_refresh_checksums.json` and
`reporting_manifest_before.csv` from immediately before the refresh.

The generated manifest differs only in the two Halistemma host-reference values.
Every publication TSV and TeX artifact is byte-identical to its pre-refresh
version. Provenance records and figure timestamps identify the new execution.
This bookkeeping correction required no sequence searches, assemblies or read
mappings. The earlier fish-rerun changes are documented separately above.

Final dry runs for Halistemma depletion, cohort remapping and expanded publication
all report nothing to do. All four affected rule-graph checks pass. Local
Markdown links and `git diff --check` pass; the three regenerated cohort PDFs
differ only in date metadata.
