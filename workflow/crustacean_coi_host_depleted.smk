"""Apply the same COI recruitment and verification to 161 host-depleted libraries.

Reuse completed recruitment/alignment for ten Nanomia libraries. Search all
retained pairs in the other 151 libraries, then verify every candidate using
the curated mitochondrial/prokaryotic panel. No original-analysis rules exist
in this bounded entry point.
"""

import csv

configfile: "config/crustacean_coi_host_depleted.json"

RUN, PANEL = config["run"], config["panel"]
PY, DIAMOND = config["python"], config["diamond"]
BLAST, BB, SEQKIT = config["blast_bin"], config["bbtools_bin"], config["seqkit"]
VERIFIED = config["verification_run"]
VERIFY_ROOT = f"{RUN}/verified"
with open(config["manifest"]) as handle:
    SAMPLES = {r["analysis_sample_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
with open(config["reuse_manifest"]) as handle:
    REUSE = {r["analysis_sample_id"] for r in csv.DictReader(handle, delimiter="\t")}

wildcard_constraints:
    sample="[^/]+",
    code="[15]",


def reads(wildcards, mate):
    return SAMPLES[wildcards.sample][f"r{mate}_path"]


def control_gate(wildcards):
    return [f"{RUN}/VERIFICATION_PASSED"]


def expected_pairs(wildcards):
    return "--expected-pairs " + SAMPLES[wildcards.sample]["pairs"]


def candidates(wildcards):
    root = config["reuse_run"] if wildcards.sample in REUSE else RUN
    return f"{root}/{wildcards.sample}"


def coi_verify_query(wildcards):
    return candidates(wildcards) + "/candidates.fna"


def coi_verify_db(wildcards):
    return f"{VERIFIED}/curated/full.dmnd"


def coi_verify_gate(wildcards):
    return [f"{RUN}/VERIFICATION_PASSED"]


include: "rules/crustacean_coi_search.smk"


rule validate_coi_expansion:
    input:
        f"{VERIFIED}/gate.json", config["manifest"],
        "config/crustacean_coi_host_depleted.json", "scripts/crustacean_coi/check_expansion.py",
    output:
        f"{RUN}/VERIFICATION_PASSED",
    shell:
        "{PY:q} scripts/crustacean_coi/check_expansion.py --config config/crustacean_coi_host_depleted.json"


include: "rules/crustacean_coi_origin_search.smk"


rule report_host_depleted_coi:
    input:
        aa=expand(f"{RUN}/verified/{{sample}}/code{{code}}.tsv", code=[1, 5], allow_missing=True),
        nt=lambda w: candidates(w) + "/competitive.tsv",
        pairs=lambda w: candidates(w) + "/candidate_pairs.tsv",
        metrics=lambda w: candidates(w) + "/recruitment.json",
        reference=f"{VERIFIED}/curated/reference_metadata.tsv",
        protein=f"{VERIFIED}/curated/protein_metadata.tsv",
        code=["scripts/crustacean_coi/verify_sample.py", "scripts/crustacean_coi/verify.py"],
        config=config["verification_config"],
    output:
        f"{RUN}/verified/{{sample}}/classified_pairs.tsv",
        f"{RUN}/verified/{{sample}}/summary.json",
    params:
        source=candidates,
        output=lambda w: f"{RUN}/verified/{w.sample}",
    resources:
        mem_mb=4000,
    shell:
        "{PY:q} scripts/crustacean_coi/verify_sample.py --directory {params.source:q} "
        "--output {params.output:q} --config {input.config:q}"


rule host_depleted:
    input:
        expand(f"{RUN}/verified/{{sample}}/summary.json", sample=SAMPLES),
    default_target: True
