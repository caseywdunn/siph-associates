"""Describe C. acnes detections alongside completed human-compatible prefix counts.

Only completed reporting tables are inputs. This entry point cannot regenerate
raw-read screens, existing mappings or catalogs.
"""

import json
import sys
from pathlib import Path

CONFIGURATION = "config/human_cutibacterium.json"
SETTINGS = json.loads(Path(CONFIGURATION).read_text())
OUT = SETTINGS["output_dir"]
SOURCES = [CONFIGURATION, "scripts/summarize_human_cutibacterium.py", "workflow/human_cutibacterium.smk"]


rule summarize_human_cutibacterium:
    input:
        sources=SOURCES,
        tables=[SETTINGS[key] for key in ("metadata", "bacterial_detections", "human_summary")],
    output:
        tables=expand(OUT + "/{name}", name=["sample_comparison.tsv", "cutibacterium_specimens.tsv", "group_summary.tsv", "summary.json", "values.tex", "provenance.json"]),
        snapshots=expand(OUT + "/provenance/source/{source}", source=SOURCES),
    params: python=sys.executable
    resources: mem_mb=1000, runtime=5
    log: "logs/human_cutibacterium.log"
    shell:
        "{params.python:q} scripts/summarize_human_cutibacterium.py --config {CONFIGURATION:q} > {log:q} 2>&1"


rule all:
    input: rules.summarize_human_cutibacterium.output
    default_target: True
