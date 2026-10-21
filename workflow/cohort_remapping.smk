"""Map 213 libraries to reviewed v3 catalogs using the shared Phase-5 rules.

Verify authorization, catalog checksums and existing read provenance; regenerate
only missing full trimmed libraries. No assembly or catalog-construction rules
are imported. The final target stops after mapping and coverage validation.
"""

import csv
import json
from pathlib import Path

configfile: "config/cohort_remapping.json"

WORK, SCRATCH = config["work_root"], config["scratch_root"]
PY = config["python"]
with open(config["mapping_read_inventory"]) as handle:
    INPUT_ROWS = list(csv.DictReader(handle, delimiter="\t"))
INPUTS = {r["sample_id"]: r for r in INPUT_ROWS}
with open(config["samples"]) as handle:
    SAMPLE_ROWS = list(csv.DictReader(handle, delimiter="\t"))
with open(config["catalog_inventory"]) as handle:
    CATALOG_FILES = [r["path"] for r in csv.DictReader(handle, delimiter="\t")]
SAMPLES = {r["sample_id"]: dict(r, read_pairs=INPUTS[r["sample_id"]]["raw_pairs"]) for r in SAMPLE_ROWS}
SAMPLE_IDS = sorted(SAMPLES)
if len(INPUT_ROWS) != len(INPUTS) or set(INPUTS) != set(SAMPLES) or len(SAMPLE_IDS) != config["expected_samples"]:
    raise WorkflowError("Mapping inventory must match all 213 unique cohort libraries")
REGENERATE = sorted(s for s in SAMPLE_IDS if INPUTS[s]["availability"] != "ready")
RETAINED_TRIMMED = sorted(s for s in SAMPLE_IDS if INPUTS[s]["availability"] == "ready" and INPUTS[s]["host_route"] == "none")
config["mapping_regenerated_samples"] = REGENERATE
PHASE5 = json.loads(Path(config["mapping_policy"]).read_text())
PHASE5["catalog_version"] = config["catalog_version"]
for resource in PHASE5["resources"].values():
    resource["partition"] = "ycga"

wildcard_constraints:
    sample="[A-Za-z0-9_.-]+",

def raw_r1(w):
    return INPUTS[w.sample]["raw_r1_paths"].split(";")

def raw_r2(w):
    return INPUTS[w.sample]["raw_r2_paths"].split(";")

def phase5_reads(w, mate):
    if INPUTS[w.sample]["availability"] == "ready":
        return INPUTS[w.sample][f"r{mate}"]
    return f"{WORK}/phase5_mapping/trimmed/{w.sample}_R{mate}.fastq.gz"

def phase5_scope_samples(w):
    if w.scope != "cohort":
        raise WorkflowError("The expanded remapping target is the full cohort")
    return SAMPLE_IDS

def phase5_resources(name):
    return {k: v for k, v in PHASE5["resources"][name].items() if k != "threads"}

def phase5_mapping_resource(key):
    def resource(w):
        size = "large" if int(INPUTS[w.sample]["read_pairs"]) >= PHASE5["large_input_pairs"] else "small"
        return PHASE5["resources"][f"mapping_{size}"][key]
    return resource

rule verify_remapping_inputs:
    input:
        config="config/cohort_remapping.json",
        authorization=config["authorization"],
        inventory=config["mapping_read_inventory"],
        samples=config["samples"],
        readiness=config["catalog_readiness"],
        catalogs=config["catalog_inventory"],
        policy=config["mapping_policy"],
        catalog_files=CATALOG_FILES,
        code=["scripts/verify_remapping_inputs.py", "scripts/assembly_completion.py"],
        provenance=[INPUTS[s]["baseline_provenance"] for s in SAMPLE_IDS],
        reads=[INPUTS[s][f"r{m}"] for s in SAMPLE_IDS if INPUTS[s]["availability"] == "ready" for m in (1, 2)],
        host_counts=[f"{SAMPLES[s]['input_root']}/phase3_cohort/inputs/{s}.json" for s in SAMPLE_IDS if INPUTS[s]["host_route"] != "none"],
    output:
        report=config["mapping_preflight"],
        eligibility=config["mapping_eligibility"],
        qc=config["mapping_library_qc"],
        retained_trim=expand(f"{WORK}/provenance/phase5_mapping/trim/{{sample}}.json", sample=RETAINED_TRIMMED),
    resources: mem_mb=4000, runtime=60, partition="ycga"
    log: f"{WORK}/logs/preflight.log"
    shell:
        "{PY:q} scripts/verify_remapping_inputs.py preflight --config {input.config:q} > {log:q} 2>&1"

include: "rules/phase5_mapping.smk"

rule verify_regenerated_mapping_reads:
    input:
        provenance=rules.trim_reads_for_mapping.output.provenance,
        preflight=rules.verify_remapping_inputs.output.report,
        baseline=lambda w: INPUTS[w.sample]["baseline_provenance"],
        r1=rules.trim_reads_for_mapping.output.r1,
        r2=rules.trim_reads_for_mapping.output.r2,
        code=["scripts/verify_remapping_inputs.py", "scripts/assembly_completion.py"],
    output: f"{WORK}/read_checks/{{sample}}.json"
    resources: mem_mb=1000, runtime=10, partition="ycga"
    log: f"{WORK}/logs/read_checks/{{sample}}.log"
    shell:
        "{PY:q} scripts/verify_remapping_inputs.py regenerated --config config/cohort_remapping.json --sample {wildcards.sample:q} > {log:q} 2>&1"

rule remapping:
    input: rules.phase5_mapping.input
    default_target: True
