"""Bounded 16S revision using the shared matching rules and completed libraries."""
import json
from pathlib import Path
ROOT = Path.cwd()
WORK = str(ROOT / "data/results")


def phase4_resources(name):
    return {"mem_mb": 16000, "runtime": 120, "partition": "ycga"}


include: "rules/cnidarian_16s.smk"

rule revision_16s:
    input: rules.cnidarian_16s.input
    default_target: True
