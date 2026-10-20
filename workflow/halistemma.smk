"""Pilot the new H. rubrum reference, deplete both libraries, then recover assemblies.

Reuse the cohort's trimming, BWA depletion and assembly/recovery rules. This
bounded workflow stops at per-library QC; it cannot build a pooled catalog or
start cohort-wide competitive mapping. Both original full-read assemblies remain.
"""

import csv
import hashlib
import json
from pathlib import Path

configfile: "config/config.yaml"
configfile: "config/halistemma.json"

REPO = Path.cwd().resolve()
WORK, SCRATCH, FIXTURE_ROOT = (config[key] for key in ("work_root", "scratch_root", "fixture_root"))
PY = config["python"]
with open(config["samples"]) as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
SAMPLES = {row["sample_id"]: row for row in rows}
SAMPLE_IDS = sorted(SAMPLES)
assert len(rows) == 2 and set(SAMPLES) == {"Ahuja2024__CWD12", "Ahuja2024__CWD17"}
assert all(row["host_route"] == "H_rubrum" for row in rows)
for resources in config["resources"].values():
    resources["partition"] = "ycga"
PHASE3_COHORT_CONFIG_PATH = REPO / "config/halistemma_cohort.json"
PHASE3_COHORT = json.loads(PHASE3_COHORT_CONFIG_PATH.read_text())
PHASE3_ASSEMBLY = json.loads((REPO / "config/phase3_assembly.json").read_text())
for resources in PHASE3_ASSEMBLY["resources"].values():
    resources["partition"] = "ycga"
PHASE3_ASSEMBLY["resources"]["assembly_large"].update(mem_mb=64000, runtime=2880)
PHASE3_ASSEMBLY_MEMBERSHIP_PATH = Path(config["assembly_gate"])
ASSEMBLY_IDS = SAMPLE_IDS
ASSEMBLY_BIN = config["tool_prefixes"]["assembly"] + "/bin"
PILOT_GATE = f"{WORK}/pilot/review.json"


def source_r1(w):
    return SAMPLES[w.sample]["r1_paths"].split(";")


def source_r2(w):
    return SAMPLES[w.sample]["r2_paths"].split(";")


raw_r1, raw_r2 = source_r1, source_r2


def host_route(w):
    return "H_rubrum"


def reference_source(w):
    assert w.route == "H_rubrum"
    return f"{WORK}/reference/GCA_987347065.1.fna"


def host_index_inputs(w):
    return [f"{WORK}/reference_indices/H_rubrum/host{suffix}" for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")] + [PILOT_GATE]


def phase3_cohort_strategy(w):
    return "reference_depleted"


def phase3_cohort_input_path(w, mate):
    return f"{WORK}/phase3_cohort/inputs/{w.sample}_R{mate}.fastq.gz"


def phase3_assembly_resource(key):
    return PHASE3_ASSEMBLY["resources"]["assembly_large"][key]


rule fetch_halistemma_reference:
    input:
        configuration="config/halistemma.json",
        report=config["host_reference"]["dataset_report"],
        inventory=config["host_reference"]["assembly_report"],
        script="scripts/halistemma/fetch_reference.py",
    output:
        reference=f"{WORK}/reference/GCA_987347065.1.fna",
        provenance=f"{WORK}/reference/validation.json",
    log: f"{WORK}/logs/reference/download.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=240,
        partition="ycga",
    shell:
        "{PY:q} {input.script:q} --config {input.configuration:q} --output {output.reference:q} --provenance {output.provenance:q} > {log:q} 2>&1"


rule freeze_halistemma_inputs:
    input:
        samples=config["samples"],
        manifest=config["manifest"],
        configuration="config/halistemma.json",
        cohort=str(PHASE3_COHORT_CONFIG_PATH),
        baselines=expand(config["baseline_trim_root"] + "/{sample}.baseline_trim.json", sample=SAMPLE_IDS),
    output: f"{WORK}/provenance/run/inputs.sha256"
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        Path(output[0]).write_text("".join(f"{hashlib.sha256(Path(p).read_bytes()).hexdigest()}  {p}\n" for p in input))


include: "rules/preprocess.smk"
include: "rules/host.smk"


rule verify_halistemma_trim:
    input:
        current=f"{WORK}/provenance/trim/{{sample}}.json",
        baseline=config["baseline_trim_root"] + "/{sample}.baseline_trim.json",
        script="scripts/halistemma/review_inputs.py",
    output: f"{WORK}/trim/{{sample}}.verified.json"
    shell:
        "{PY:q} {input.script:q} trim --current {input.current:q} --baseline {input.baseline:q} --output {output:q}"


rule sample_halistemma_pilot:
    input:
        r1=f"{SCRATCH}/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{SCRATCH}/trimmed/{{sample}}_R2.fastq.gz",
        validated=f"{WORK}/trim/{{sample}}.verified.json",
    output:
        r1=f"{WORK}/pilot/{{sample}}_R1.fastq.gz",
        r2=f"{WORK}/pilot/{{sample}}_R2.fastq.gz",
    log: f"{WORK}/logs/pilot/{{sample}}.sampling.log"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=240,
        partition="ycga",
    conda: "../envs/fastp.yaml"
    shell:
        """
        reformat.sh -Xmx12g in1={input.r1:q} in2={input.r2:q} \
          out1={output.r1:q} out2={output.r2:q} \
          samplereadstarget={config[pilot_pairs]} sampleseed={config[pilot_seed]} \
          threads={threads} overwrite=t > {log:q} 2>&1
        """


rule map_halistemma_pilot:
    input:
        r1=f"{WORK}/pilot/{{sample}}_R1.fastq.gz",
        r2=f"{WORK}/pilot/{{sample}}_R2.fastq.gz",
        index=expand(f"{WORK}/reference_indices/H_rubrum/host{{suffix}}", suffix=[".amb", ".ann", ".bwt", ".pac", ".sa"]),
    output: f"{WORK}/pilot/{{sample}}.bam"
    params:
        bwa=f"{ASSEMBLY_BIN}/bwa",
        samtools=f"{ASSEMBLY_BIN}/samtools",
        index=f"{WORK}/reference_indices/H_rubrum/host",
    log: f"{WORK}/logs/pilot/{{sample}}.mapping.log"
    benchmark: f"{WORK}/benchmarks/pilot/{{sample}}.tsv"
    threads: 16
    resources:
        mem_mb=64000,
        runtime=120,
        partition="ycga",
    shell:
        """
        {params.bwa:q} mem -t {threads} {params.index:q} {input.r1:q} {input.r2:q} 2> {log:q} \
          | {params.samtools:q} view -b -o {output:q}.tmp - 2>> {log:q}
        mv {output:q}.tmp {output:q}
        """


rule summarize_halistemma_pilot:
    input:
        r1=f"{WORK}/pilot/{{sample}}_R1.fastq.gz",
        r2=f"{WORK}/pilot/{{sample}}_R2.fastq.gz",
        bam=f"{WORK}/pilot/{{sample}}.bam",
        script="scripts/halistemma/review_inputs.py",
    output: f"{WORK}/pilot/{{sample}}.json"
    params: samtools=f"{ASSEMBLY_BIN}/samtools"
    resources:
        mem_mb=2000,
        runtime=60,
        partition="ycga",
    shell:
        "{PY:q} {input.script:q} pilot --r1 {input.r1:q} --r2 {input.r2:q} --bam {input.bam:q} --samtools {params.samtools:q} --sample {wildcards.sample:q} --pairs {config[pilot_pairs]} --output {output:q}"


rule review_halistemma_pilots:
    input:
        reports=expand(f"{WORK}/pilot/{{sample}}.json", sample=SAMPLE_IDS),
        samples=config["samples"],
        script="scripts/halistemma/review_inputs.py",
    output: PILOT_GATE
    shell:
        "{PY:q} {input.script:q} gate --reports {input.reports:q} --samples {input.samples:q} --minimum {config[pilot_minimum_host_pair_fraction]} --output {output:q}"


include: "rules/phase3_cohort.smk"
include: "rules/phase3_assembly.smk"


rule review_halistemma_assembly_inputs:
    input:
        eligibility=f"{WORK}/phase3_cohort/input_eligibility.tsv",
        samples=config["samples"],
        script="scripts/halistemma/review_inputs.py",
    output: str(PHASE3_ASSEMBLY_MEMBERSHIP_PATH)
    params: minimum=PHASE3_COHORT["minimum_assembly_input_pairs"]
    shell:
        "{PY:q} {input.script:q} membership --eligibility {input.eligibility:q} --samples {input.samples:q} --minimum {params.minimum} --output {output:q}"


rule halistemma_pilot:
    input: PILOT_GATE


rule halistemma:
    input: f"{WORK}/stages/phase3_assembly.done"
    default_target: True
