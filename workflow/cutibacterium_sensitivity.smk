"""Check the two primary community results after omitting one skin reference.

Completed tables are the only biological inputs. Original results stay unchanged.
"""

import json
import sys
from pathlib import Path

CONFIGURATION = "config/cutibacterium_sensitivity.json"
SETTINGS = json.loads(Path(CONFIGURATION).read_text())
OUT = SETTINGS["output_dir"]
SOURCES = [CONFIGURATION, SETTINGS["original_config"],
           "workflow/scripts/analyze_coauthor_communities.R",
           "workflow/scripts/analyze_cutibacterium_sensitivity.R",
           "scripts/run_cutibacterium_sensitivity.py",
           "workflow/cutibacterium_sensitivity.smk", "envs/stats.yaml"]


rule cutibacterium_community_sensitivity:
    input:
        sources=SOURCES,
        tables=[SETTINGS[key] for key in ("metadata", "membership", "original_tests", "bacterial_grades")],
    output:
        tables=expand(OUT + "/{name}", name=["community_tests.tsv", "community_membership.tsv", "summary.json", "values.tex", "session_info.txt", "provenance.json"]),
        snapshots=expand(OUT + "/provenance/source/{source}", source=SOURCES),
    params: python=sys.executable
    threads: 1
    resources: mem_mb=4000, runtime=10
    log: "logs/cutibacterium_sensitivity.log"
    shell:
        "{params.python:q} scripts/run_cutibacterium_sensitivity.py --config {CONFIGURATION:q} > {log:q} 2>&1"


rule all:
    input: rules.cutibacterium_community_sensitivity.output
    default_target: True
