"""Describe sampling, test identifiable communities, then record the result.

1. Join reviewed tissue/depth metadata to existing detections: prepare_coauthor_metadata.
2. Test same-tissue/host/region contrasts and plot depth: analyze_coauthor_communities.
3. Verify and preserve input/output/source identities: record_coauthor_metadata.
"""

import json
import sys
from pathlib import Path

COAUTHOR_METADATA_CONFIG = "config/coauthor_metadata_analysis.json"
COAUTHOR_METADATA = json.loads(Path(COAUTHOR_METADATA_CONFIG).read_text())
COAUTHOR_OUT = COAUTHOR_METADATA["output_dir"]
COAUTHOR_INPUTS = [
    "figures/cohort/sample_metadata.tsv", "figures/cohort/bacterial_detections.tsv",
    "figures/cohort/viral_detections.tsv", "figures/cohort/eukaryote_evidence.tsv",
    "figures/cohort/library_flowcells.tsv", "figures/cohort/mycoplasmatales_incidence.tsv",
    COAUTHOR_METADATA["mapping_summary"], COAUTHOR_METADATA["bacterial_grades"],
    COAUTHOR_METADATA["fish_summary"],
    COAUTHOR_METADATA["host_depletion_stats"], COAUTHOR_METADATA["host_depletion_provenance"],
    COAUTHOR_METADATA["historical_host_concordance"],
]
COAUTHOR_SOURCES = [
    COAUTHOR_METADATA_CONFIG, "scripts/analyze_coauthor_metadata.py",
    "scripts/analyze_collection_context.py", "workflow/scripts/analyze_coauthor_communities.R",
    "scripts/record_coauthor_metadata.py", "workflow/coauthor_metadata.smk",
    "workflow/rules/coauthor_metadata.smk", "envs/stats.yaml",
]
COAUTHOR_TABLES = [
    "specimen_context.tsv", "associate_detections.tsv", "tissue_incidence.tsv",
    "depth_incidence.tsv", "host_depth_design.tsv", "host_locality_incidence.tsv",
    "targeted_incidence.tsv", "contrast_design.tsv", "contrast_membership.tsv",
    "descriptive_values.json",
    "host_depletion_summary.tsv", "historical_host_concordance.tsv", "host_reference_values.json",
]
COAUTHOR_STATS = [
    "community_tests.tsv", "dispersion_tests.tsv", "community_membership.tsv",
    "community_values.json", "depth_patterns.pdf", "depth_patterns.png",
    "targeted_incidence.pdf", "targeted_incidence.png", "session_info.txt",
]


rule prepare_coauthor_metadata:
    input:
        code=["scripts/analyze_coauthor_metadata.py", "scripts/analyze_collection_context.py"],
        config=COAUTHOR_METADATA_CONFIG,
        tables=COAUTHOR_INPUTS,
    output:
        expand(COAUTHOR_OUT + "/{name}", name=COAUTHOR_TABLES),
    resources: mem_mb=2000, runtime=10, partition="ycga"
    params: python=sys.executable
    log: "logs/coauthor_metadata/prepare.log"
    shell:
        "{params.python:q} scripts/analyze_coauthor_metadata.py --config {input.config:q} > {log:q} 2>&1"


rule analyze_coauthor_communities:
    input:
        code="workflow/scripts/analyze_coauthor_communities.R",
        config=COAUTHOR_METADATA_CONFIG,
        tables=rules.prepare_coauthor_metadata.output,
        grades=COAUTHOR_METADATA["bacterial_grades"],
    output:
        expand(COAUTHOR_OUT + "/{name}", name=COAUTHOR_STATS),
    resources: mem_mb=4000, runtime=30, partition="ycga"
    params: rscript="/gpfs/gibbs/project/dunn/cwd7/conda_envs/stats_env/bin/Rscript"
    threads: 1
    log: "logs/coauthor_metadata/statistics.log"
    shell:
        "{params.rscript:q} {input.code:q} {input.config:q} > {log:q} 2>&1"


rule record_coauthor_metadata:
    input:
        sources=COAUTHOR_SOURCES,
        inputs=COAUTHOR_INPUTS,
        tables=rules.prepare_coauthor_metadata.output,
        statistics=rules.analyze_coauthor_communities.output,
    output:
        provenance=COAUTHOR_OUT + "/provenance.json",
        values=COAUTHOR_OUT + "/values.tex",
        snapshot=expand(COAUTHOR_OUT + "/provenance/source/{path}", path=COAUTHOR_SOURCES),
    params: python=sys.executable
    resources: mem_mb=2000, runtime=5, partition="ycga"
    shell:
        "{params.python:q} scripts/record_coauthor_metadata.py --config {COAUTHOR_METADATA_CONFIG:q}"
