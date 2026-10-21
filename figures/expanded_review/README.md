# Expanded-cohort comparison

This review compares the publication tables at revision `44ade17` with the
213-library v3 analysis. **The SSU and publication exports are complete and
validated.** Execution and reuse boundaries are in the
[reporting guide](../../docs/expanded_reporting.md).

The comparison producer translates old catalog identifiers through the frozen
v3 crosswalk, including merged viral clusters. A detection is a library/target
pair (or a library/reporting-group pair for eukaryotes), not an organism count.

## Findings established from the completed mapping and statistical tables

| Quantity | Previous report | Expanded v3 |
| --- | ---: | ---: |
| Libraries | 205 | 213 |
| Retained bacterial MAGs | 103 | 104 |
| MAG clusters | 31 | 32 |
| Bacterial detections | 557 | 567 |
| Archaeal detections | 0 | 1 |
| Libraries with a microbial detection | 138 | 140 |
| Microbial detections with assembly support | 302 | 303 |
| Associate-virus detections | 530 | 631 |
| Libraries with an associate-virus detection | 122 | 127 |
| Associate-virus detections with assembly support | 348 | 442 |
| Nonredundant eukaryotic detections | 61 | 64 |
| Libraries with a eukaryotic detection | 39 | 40 |
| Eukaryotic detections with assembly support and passing reads | 42 | 42 |
| Strong fish 12S paired sequences | 32 | 32 |
| Fish-positive libraries | 8 | 8 |

All original microbial detections remain. Nine added calls occur in Bargmannia,
all supported by reads without passing assembly support. The other two added
calls are the new Tyloplasma-related MAGSP0033/BAC00459 in Physalia
YPM-IZ-111016 and YPM-IZ-111017. Four original microbial evidence grades change;
the two reductions for MAGSP0006 follow removal of host/eukaryotic fragments
during catalog curation. Both retain read-supported presence.

The archaeal call is a Nitrosopelagicus reference match in Bargmannia stenotes
D666-D10 (GCA_013390375.1/BAC00450; 13.1901% breadth, 1,083 retained reads).
It has no passing assembly support and is below the 20% breadth sensitivity
cutoff. No archaeal MAG is recovered. The rejected MAGSP0026 remains excluded.

Viral changes comprise 13 lost calls on six removed targets, 15 gained calls in
original libraries, and 99 calls in the added libraries. Of the latter, 86 occur
in formalin-preserved Bargmannia elongata T676-SS3. Retain this specimen context
when interpreting the increase. The new fish screens add 1,430,751,915 filtered
pairs but no strong fish matches. With the two Halistemma entries selected from
host-depleted reruns, the combined fish survey contains 9,313,145,671 pairs.

T676-SS3 also has 288 insect and 71 plant SSU pairs at ≥97% identity. It is
the only library exceeding the locked 37-pair threshold in these predefined
background groups. The primary calls remain unchanged; the comparison includes
a specimen flag and a descriptive sensitivity summary omitting this library.
The background signal qualifies interpretation of its 86 viral detections,
without establishing that they are false. The paper describes a distinct DNA
extraction protocol for this formalin-preserved nectophore; the source of the
background signal is unresolved.

Omitting T676-SS3 gives 567 microbial detections in 139 libraries, 545 viral
detections in 126 libraries, and 61 eukaryotic detections in 39 libraries.
All 86 viral targets detected in T676-SS3 are absent from other libraries at
the accepted threshold. Its three added eukaryotic groups are fungi,
Bacillariophytina and the perkinsid group A31, all read-only calls. Original
eukaryotic detections and their grades are unchanged.

The seven focal Metamycoplasmataceae genomes retain identical source and sequence
hashes. Their phylogeny, functional annotations and external cnidarian comparisons
remain valid within that scope. Tyloplasma is outside those seven-genome analyses.
The eight Waki-comparison source assemblies are also unchanged.

## Outputs and handoff

`compare_expanded_publication` produces [summary.json](summary.json),
[detection_changes.tsv](detection_changes.tsv),
[original_metadata_changes.tsv](original_metadata_changes.tsv),
[retained_focal_genomes.tsv](retained_focal_genomes.tsv),
[ssu_background_flags.tsv](ssu_background_flags.tsv) and
[background_sensitivity.tsv](background_sensitivity.tsv). The original
205 libraries' reviewed metadata are unchanged. The preserved baseline is
`data/results/post_expansion/2026-10-21/reporting/baseline/figures-44ade17.tar.gz`;
the generated summary records its SHA-256.

The [manuscript handoff](../../dev_docs/expanded-manuscript-handoff.md) identifies
affected text and asset-source paths. No manuscript files were edited in this
reporting refresh. Historical v2 read-handling controls keep their original
catalog and subsets; they are not relabelled as v3 tests.
