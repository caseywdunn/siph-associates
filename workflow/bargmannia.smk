"""Process the eight newly public Manko2026 libraries with the cohort methods."""

import csv
import hashlib
import json
from pathlib import Path

configfile: "config/config.yaml"
configfile: "config/bargmannia.json"

REPO = Path.cwd().resolve()
WORK = config["work_root"]
SCRATCH = config["scratch_root"]
FIXTURE_ROOT = config["fixture_root"]
with open(config["samples"]) as handle:
    sample_rows = list(csv.DictReader(handle, delimiter="\t"))
SAMPLES = {r["analysis_sample_id"]: r for r in sample_rows}
SAMPLE_IDS = sorted(SAMPLES)
assert len(SAMPLES) == len(sample_rows) == 8
assert {r["study"] for r in sample_rows} == {"Manko2026"}
assert all(r["host_route"] == "none" for r in sample_rows)
for resources in config["resources"].values():
    resources["partition"] = "ycga"
PHASE3_COHORT_CONFIG_PATH = REPO / "config/phase3_cohort.json"
PHASE3_COHORT = json.loads(PHASE3_COHORT_CONFIG_PATH.read_text())
PHASE3_ASSEMBLY = json.loads((REPO / "config/phase3_assembly.json").read_text())
for settings in [PHASE3_COHORT, PHASE3_ASSEMBLY]:
    for resources in settings["resources"].values():
        resources["partition"] = "ycga"
PHASE3_ASSEMBLY_MEMBERSHIP_PATH = Path(config["assembly_gate"])
ASSEMBLY_IDS = SAMPLE_IDS


def source_r1(w):
    return SAMPLES[w.sample]["r1_paths"].split(";")


def source_r2(w):
    return SAMPLES[w.sample]["r2_paths"].split(";")


raw_r1, raw_r2 = source_r1, source_r2


def host_route(w):
    return "none"


def host_index_inputs(w):
    return []


def phase3_cohort_strategy(w):
    return "fixed_effort_trimmed"


def phase3_cohort_input_path(w, mate):
    return f"{WORK}/phase3_cohort/inputs/{w.sample}_R{mate}.fastq.gz"


def phase3_assembly_resource(key):
    return PHASE3_ASSEMBLY["resources"]["assembly_large"][key]


rule verify_bargmannia_raw_files:
    input:
        files=config["raw_files"],
        r1=source_r1,
        r2=source_r2,
        script="scripts/bargmannia/verify_raw.py",
    output:
        f"{WORK}/verified/{{sample}}.json",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=120,
        partition="ycga",
    shell:
        "{config[python]:q} {input.script:q} --files {input.files:q} --sample {wildcards.sample:q} --output {output:q}"


rule freeze_bargmannia_inputs:
    input:
        samples=config["samples"],
        manifest=config["manifest"],
        configuration="config/bargmannia.json",
        verified=expand(f"{WORK}/verified/{{sample}}.json", sample=SAMPLE_IDS),
    output:
        f"{WORK}/provenance/run/inputs.sha256",
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        Path(output[0]).write_text("".join(f"{hashlib.sha256(Path(p).read_bytes()).hexdigest()}  {p}\n" for p in input))


include: "rules/preprocess.smk"
include: "rules/screens.smk"
include: "rules/phase3_cohort.smk"
include: "rules/phase3_assembly.smk"


rule review_bargmannia_assembly_inputs:
    input:
        f"{WORK}/phase3_cohort/input_eligibility.tsv",
    output:
        str(PHASE3_ASSEMBLY_MEMBERSHIP_PATH),
    run:
        with open(input[0]) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert {r["sample_id"] for r in rows} == set(SAMPLE_IDS)
        if any(r["assembly_eligible"] != "true" for r in rows):
            raise ValueError("A new library falls below the established assembly threshold; review before proceeding")
        fields = ["sample_id", "host_route", "strategy", "retained_pairs", "minimum_assembly_input_pairs", "assembly_eligible"]
        with open(output[0], "w") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow(dict(row, minimum_assembly_input_pairs=PHASE3_COHORT["minimum_assembly_input_pairs"]))


# The same fish filter/alignment rules consume the retained full trimmed reads.
ROOT = config["run"]
PY = config["python"]
READ_LABEL = config["read_pair_id"]


def fish_read(w, mate):
    return f"{SCRATCH}/trimmed/{w.sample}_R{mate}.fastq.gz"


def fish_manifest(w):
    return f"{ROOT}/inputs/{w.sample}.tsv"


def fish_index(w):
    return 0


rule inventory_bargmannia_fish_reads:
    input:
        r1=lambda w: fish_read(w, 1),
        r2=lambda w: fish_read(w, 2),
        provenance=f"{WORK}/provenance/trim/{{sample}}.json",
        report=f"{WORK}/trim/{{sample}}.fastp.json",
        samples=config["samples"],
        script="scripts/bargmannia/read_inventory.py",
    output:
        f"{ROOT}/inputs/{{sample}}.tsv",
    shell:
        "{PY:q} {input.script:q} --samples {input.samples:q} --sample {wildcards.sample:q} --root {WORK:q} --output {output:q}"


include: "rules/fish_12s_search.smk"


rule review_bargmannia_fish:
    input:
        hits=f"{ROOT}/union/{{sample}}/{READ_LABEL}/alignment/hits.tsv",
        manifest=fish_manifest,
        metadata=config["reference_metadata"],
    output:
        f"{ROOT}/review/{{sample}}/sample_summary.tsv",
    params:
        output=lambda w: f"{ROOT}/review/{w.sample}",
    shell:
        "{PY:q} scripts/fish_12s/summarize_cohort.py --config config/bargmannia.json --manifest {input.manifest:q} --reference-metadata {input.metadata:q} --output {params.output:q}"


rule bargmannia:
    input:
        f"{WORK}/stages/phase3_assembly.done",
        expand(f"{WORK}/screens/bracken/{{sample}}.S.tsv", sample=SAMPLE_IDS),
        expand(f"{WORK}/screens/sylph/{{sample}}.profile.tsv", sample=SAMPLE_IDS),
        expand(f"{WORK}/screens/phyloflash/{{sample}}.tar.gz", sample=SAMPLE_IDS),
        expand(f"{ROOT}/review/{{sample}}/sample_summary.tsv", sample=SAMPLE_IDS),
    default_target: True
