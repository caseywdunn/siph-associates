"""Extend the tested fish filters to the remaining 156 host-depleted libraries.

Use the pilot implementations unchanged, retaining both mates of their union.
The five completed fish pilots remain external inputs, without repeating scans.
Every new union receives the same competitive whole-mitochondrial BLAST search.
"""

import csv

configfile: "config/fish_12s_host_depleted_cohort.json"

ROOT, PY = config["run"], config["python"]
with open(config["manifest"]) as handle:
    MANIFEST = list(csv.DictReader(handle, delimiter="\t"))
SAMPLES = {r["analysis_sample_id"]: dict(r, index=i) for i, r in enumerate(MANIFEST)}
REUSE = set()
if config.get("reuse_manifest"):
    with open(config["reuse_manifest"]) as handle:
        REUSE = {r["analysis_sample_id"] for r in csv.DictReader(handle, delimiter="\t")}
NEW = sorted(SAMPLES.keys() - REUSE)

wildcard_constraints:
    sample="[^/]+",


READ_LABEL = "host_depleted"


def fish_read(w, mate):
    return SAMPLES[w.sample][f"r{mate}_path"]


def fish_manifest(w):
    return config["manifest"]


def fish_index(w):
    return SAMPLES[w.sample]["index"]


include: "rules/fish_12s_search.smk"


rule validate_fish_cohort_controls:
    input:
        f"{ROOT}/union/synthetic/host_depleted/alignment/hits.tsv",
        f"{ROOT}/union/synthetic/host_depleted/candidates.fasta",
        f"{ROOT}/union/synthetic/host_depleted/candidate_pairs.tsv",
        f"{ROOT}/union/synthetic/host_depleted/run.json",
        truth="data/results/fish_12s/positive_aligner_controls/2026-10-17/synthetic/known_answer/checks.json",
        code="scripts/fish_12s/check_cohort_controls.py",
        classifier="scripts/fish_12s/classify_candidates.py",
        metadata="data/results/fish_12s/reference_metadata/2026-10-18/references.tsv",
    output:
        f"{ROOT}/CONTROLS_PASSED",
    shell:
        "{PY:q} {input.code:q} --root {ROOT:q} --truth {input.truth:q} --reference-metadata {input.metadata:q}"


rule host_depleted:
    input:
        expand(f"{ROOT}/union/{{sample}}/host_depleted/alignment/hits.tsv", sample=NEW),
        expand(config["reuse_screen"] + "/{sample}/host_depleted/alignment/hits.tsv", sample=REUSE),
        expand(config["reuse_aligners"] + "/{sample}/host_depleted/run.json", sample=REUSE),
        ["data/results/fish_12s/host_depleted_aligner_only/2026-10-17/Church2025__YPM-IZ-110907/host_depleted/alignment/hits.tsv"] if REUSE else [],
        [f"{ROOT}/CONTROLS_PASSED"] if config.get("synthetic_validation") else [],
    default_target: True
