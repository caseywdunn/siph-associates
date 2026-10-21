# Halistemma host-depleted fish screens, 2026-10-21

Casey Dunn requested that both Halistemma rubrum fish screens use all retained
host-depleted pairs, matching Physalia and Nanomia septata, and that the root
README give fish screening its own row in the analysis-steps table. Work is on
`main`, from source checkpoint `043a2bb` with the source changes captured in
the execution archives. No microbial/viral analysis or host mapping is repeated.

## Inputs and execution

The two input rows in `config/fish_12s_halistemma_host_depleted.tsv` were derived
from `data/results/halistemma/2026-10-19/phase3_cohort/inputs/{sample}.json`.
File sizes were checked against the corresponding existing depletion provenance.
Both use `H_rubrum` (GCA_987347065.1) and require both mates unmapped; there is
no additional cap or subsampling. Original raw-read assignments are unchanged.

| Library | Pairs before host depletion | All retained pairs screened |
| --- | ---: | ---: |
| Ahuja2024__CWD12 | 100,871,052 | 131,703 |
| Ahuja2024__CWD17 | 165,797,762 | 2,352,421 |

Allocation `13691241` used the existing host-depleted entry point with
`config/fish_12s_halistemma_host_depleted.json` and target
`review_host_depleted_fish`: 24 CPUs, 72 GB, four-hour limit in YCGA. All nine
steps completed successfully. It reused the original masked fish reference,
Bowtie2 index, competitive mitochondrial database and accession-linked taxonomy;
the launcher's reference checks and fish/non-fish control gate passed. Regex,
BBDuk, Bowtie2, SortMeRNA, paired union recovery and competitive BLAST commands
and scientific thresholds are unchanged. The union step checks synchronization
and exact input-pair accounting.

Outputs: `data/results/fish_12s/halistemma_host_depleted/2026-10-21/`.
Execution logs, tool versions, source archives and checksums are under
`data/results/fish_12s/executions/13691241/`. The prior fish publication bundle
is preserved as `baseline/fish_report_043a2bb.tar.gz` beneath the output root;
the original full-read searches and reports remain in place.

Dependent allocation `13691335` completed the reporting-only
`expanded_fish_summary` target (two CPUs, 8 GB, 30-minute limit). Its execution
record is under `data/results/post_expansion/2026-10-21/reporting/executions/13691335/`.
Selection excludes the two old full-read reports, includes their new depleted
report, and checks 213 unique library IDs. No other fish search was repeated.

## Results and checks

CWD12 has no recruited candidates. CWD17 has three candidate pairs, all in the
`other_or_no_hit` category. Neither has a qualifying or weaker fish-similar pair.
The combined report contains 163 depleted and 50 full-read libraries, totaling
9,313,145,671 pairs (419,225,871 depleted; 8,893,919,800 full-read).
The 32 distinct qualifying pairs in eight Physalia libraries, their recruitment
method counts, and the 81 weaker fish-similar pairs are unchanged.

The 13 classification/report-selection tests passed. The four reference-assessment
tests also passed in their documented Arrow/Python environment; their initial
attempt in the Snakemake environment lacked PyArrow. No reference or environment
was changed. The test for report selection rejects duplicate and unmatched
replacement libraries. Current launch instructions and the generated graph are
in the [fish workflow guide](../../docs/prey_marker_screen.md).

Direct comparison with the archived publication bundle confirms that the other
211 sample-summary rows are unchanged and `positive_pairs.tsv` is byte-identical.
The baseline archive SHA-256 is
`d4cc649ec76b2b5756c964805125411028425df7fc62fbded3677659f73420c4`.
Completed-run dry runs for `review_host_depleted_fish` and
`expanded_publication` both report that all requested files are up to date.
Local Markdown links and whitespace checks pass. The affected rule graphs were
regenerated from the workflows.
