"""Load the existing catalog rules for bounded reuse during finalization."""

import json
from pathlib import Path

WORK, SCRATCH = config["work_root"], config["scratch_root"]
PHASE3_CATALOG_CONFIG_PATH = Path("config/phase3_catalog.json")
PHASE4_CATALOG_CONFIG_PATH = Path("config/catalog_freeze_policy.json")
PHASE3_CATALOG = json.loads(PHASE3_CATALOG_CONFIG_PATH.read_text())
PHASE4_CATALOG = json.loads(PHASE4_CATALOG_CONFIG_PATH.read_text())
PHASE3_ASSEMBLY_MEMBERSHIP_PATH = Path(config["preparation_root"]) / "phase3_cohort/membership.tsv"
for settings in (PHASE3_CATALOG, PHASE4_CATALOG):
    for resources in settings["resources"].values():
        resources["partition"] = "ycga"

def catalog_resources(name):
    return {k: v for k, v in PHASE3_CATALOG["resources"][name].items() if k != "threads"}

def phase4_resources(name):
    return {k: v for k, v in PHASE4_CATALOG["resources"][name].items() if k != "threads"}

include: "rules/phase3_catalog.smk"
include: "rules/phase4_catalog.smk"
