"""Finish per-library assembly and prepare QC inputs; no cohort-mapping rules."""

import csv
import json
from pathlib import Path

configfile: "config/assembly_completion.json"

WORK = config["work_root"]
# The shared binning rule creates its archive from inside the bins directory.
# Its temporary destination must remain absolute after that directory change.
SCRATCH = str((Path(WORK) / "scratch").resolve())
PY = config["python"]
PHASE3_ASSEMBLY_MEMBERSHIP_PATH = Path(config["membership"])
with PHASE3_ASSEMBLY_MEMBERSHIP_PATH.open() as handle:
    COMPLETION_ROWS = list(csv.DictReader(handle, delimiter="\t"))
COMPLETION = {r["sample_id"]: r for r in COMPLETION_ROWS}
ASSEMBLY_IDS = sorted(COMPLETION)
assert len(COMPLETION_ROWS) == len(COMPLETION) == config["expected_completion_samples"]
assert all(r["assembly_eligible"] == "true" and int(r["minimum_assembly_input_pairs"]) == 0 for r in COMPLETION_ROWS)
PHASE3_ASSEMBLY = json.loads(Path("config/phase3_assembly.json").read_text())
for resources in PHASE3_ASSEMBLY["resources"].values():
    resources["partition"] = "ycga"
config["assembly_gate"] = f"{WORK}/phase3_cohort/input_eligibility.tsv"


def phase3_cohort_input_path(w, mate):
    return f"{COMPLETION[w.sample]['input_root']}/phase3_cohort/inputs/{w.sample}_R{mate}.fastq.gz"


def phase3_assembly_resource(key):
    return PHASE3_ASSEMBLY["resources"]["assembly_small"][key]


rule validate_completion_inputs:
    input:
        settings="config/assembly_completion.json",
        membership=str(PHASE3_ASSEMBLY_MEMBERSHIP_PATH),
        code="scripts/assembly_completion.py",
        reads=[f"{r['input_root']}/phase3_cohort/inputs/{r['sample_id']}_R{m}.fastq.gz" for r in COMPLETION_ROWS for m in (1, 2)],
        records=[f"{r['input_root']}/{part}/{r['sample_id']}.json" for r in COMPLETION_ROWS for part in ("phase3_cohort/inputs", "provenance/phase3_cohort/inputs", "provenance/trim")],
    output:
        eligibility=f"{WORK}/phase3_cohort/input_eligibility.tsv",
        provenance=f"{WORK}/provenance/selected_inputs.json",
    resources:
        mem_mb=1000,
        runtime=10,
        partition="ycga",
    shell:
        "{PY:q} {input.code:q} validate --config {input.settings:q}"


include: "rules/phase3_assembly.smk"


rule assemblies:
    input: f"{WORK}/stages/phase3_assembly.done"
    default_target: True
