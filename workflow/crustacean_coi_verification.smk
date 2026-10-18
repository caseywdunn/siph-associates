"""Audit COI reference origins, then recheck preserved reads without library scans.

Annotated mitochondrial COI and UniProt oxidases anchor reference-label review.
Curated protein panels verify real candidates and original/held-out controls.
Nucleotide evidence is reused from the completed pilot and sensitivity benchmark.
"""

configfile: "config/crustacean_coi_verification.json"

from pathlib import Path

ROOT = config["run"]
SOURCE = config["source_panel"]
PY, DIAMOND = config["python"], config["diamond"]
CONFIG = "config/crustacean_coi_verification.json"
CODE = "scripts/crustacean_coi/verify.py"
PILOT = config["pilot_run"]
VALIDATION = config["validation_run"]
PILOT_SAMPLES = sorted(p.name for p in Path(PILOT).glob("Ahuja2026__*"))
DATASETS = ["pilot", "controls", "heldout"]
VERIFY_ROOT = f"{ROOT}/verification"


def coi_verify_query(wildcards):
    return f"{ROOT}/queries/{wildcards.sample}.fna"


def coi_verify_db(wildcards):
    panel = "heldout" if wildcards.sample == "heldout" else "full"
    return f"{ROOT}/curated/{panel}.dmnd"


def coi_verify_gate(wildcards):
    return []

wildcard_constraints:
    dataset="pilot|controls|heldout",
    panel="full|heldout",
    code="[15]",


rule fetch_prokaryotic_oxidases:
    input:
        "scripts/crustacean_coi/fetch_oxidases.py", CONFIG,
    output:
        f"{ROOT}/source/uniprot_oxidases.tsv", f"{ROOT}/source/uniprot_receipt.json",
    shell:
        "{PY:q} scripts/crustacean_coi/fetch_oxidases.py --config {CONFIG:q}"


rule prepare_oxidase_anchors:
    input:
        rules.fetch_prokaryotic_oxidases.output,
        config["mitochondrial_genbank"], f"{SOURCE}/reference_metadata.tsv",
        f"{SOURCE}/annotated.faa", CODE, CONFIG,
    output:
        f"{ROOT}/anchors/reference.faa", f"{ROOT}/anchors/reference_metadata.tsv",
        f"{ROOT}/anchors/protein_metadata.tsv", f"{ROOT}/anchors/summary.json",
    resources:
        mem_mb=3000,
    shell:
        "{PY:q} {CODE:q} anchors --config {CONFIG:q}"


rule index_oxidase_anchors_diamond:
    input:
        f"{ROOT}/anchors/reference.faa",
    output:
        f"{ROOT}/anchors/reference.dmnd",
    threads: 4
    resources:
        mem_mb=2000,
    shell:
        "{DIAMOND:q} makedb --in {input:q} --db {output:q} --threads {threads}"


rule audit_target_reference_origins_diamond:
    input:
        query=f"{ROOT}/anchors/target_unique.faa", db=f"{ROOT}/anchors/reference.dmnd",
    output:
        f"{ROOT}/anchors/target_audit.tsv",
    log:
        f"{ROOT}/anchors/target_audit.log",
    threads: 8
    resources:
        mem_mb=6000,
    benchmark:
        f"{ROOT}/anchors/benchmark.tsv",
    shell:
        "{DIAMOND:q} blastp --query {input.query:q} --db {input.db:q} --out {output:q} "
        "--very-sensitive --evalue 1e-5 --query-cover 70 --top 5 --max-hsps 1 "
        "--block-size 0.4 --index-chunks 1 --threads {threads} "
        "--outfmt 6 qseqid sseqid pident length qcovhsp bitscore evalue 2> {log:q}"


rule curate_coi_assignment_references:
    input:
        rules.prepare_oxidase_anchors.output,
        rules.audit_target_reference_origins_diamond.output,
        f"{SOURCE}/target.faa", f"{SOURCE}/target.fna",
        f"{ROOT}/anchors/target_representatives.tsv",
        f"{VALIDATION}/panel/reference_metadata.tsv", config["manual_quarantine"], CODE, CONFIG,
    output:
        expand(f"{ROOT}/curated/{{name}}", name=["reference_audit.tsv", "reference_metadata.tsv", "protein_metadata.tsv", "full.faa", "heldout.faa", "quarantined.fna", "summary.json"]),
    resources:
        mem_mb=4000,
    shell:
        "{PY:q} {CODE:q} curate --config {CONFIG:q}"


rule deduplicate_target_audit_queries:
    input:
        f"{SOURCE}/target.faa", CODE, CONFIG,
    output:
        f"{ROOT}/anchors/target_unique.faa", f"{ROOT}/anchors/target_representatives.tsv",
    resources:
        mem_mb=2000,
    shell:
        "{PY:q} {CODE:q} deduplicate --config {CONFIG:q}"


rule index_curated_coi_diamond:
    input:
        f"{ROOT}/curated/{{panel}}.faa",
    output:
        f"{ROOT}/curated/{{panel}}.dmnd",
    threads: 4
    resources:
        mem_mb=2000,
    shell:
        "{DIAMOND:q} makedb --in {input:q} --db {output:q} --threads {threads}"


rule prepare_preserved_coi_candidates:
    input:
        expand(f"{PILOT}/{{sample}}/{{name}}", sample=PILOT_SAMPLES + ["controls"], name=["candidates.fna", "candidate_pairs.tsv", "competitive.tsv"]),
        expand(f"{VALIDATION}/heldout/{{name}}", name=["candidates.fna", "candidate_pairs.tsv", "competitive.tsv"]),
        CODE, CONFIG,
    output:
        expand(f"{ROOT}/queries/{{dataset}}.{{suffix}}", dataset=DATASETS, suffix=["fna", "pairs.tsv", "nt.tsv"]),
    shell:
        "{PY:q} {CODE:q} queries --config {CONFIG:q}"


include: "rules/crustacean_coi_origin_search.smk"


rule report_coi_origin_verification:
    input:
        expand(f"{ROOT}/verification/{{dataset}}/code{{code}}.tsv", dataset=DATASETS, code=[1, 5]),
        rules.prepare_preserved_coi_candidates.output,
        rules.curate_coi_assignment_references.output,
        f"{PILOT}/controls/truth.json", f"{VALIDATION}/heldout/truth.json", CODE, CONFIG,
    output:
        f"{ROOT}/summary.json", f"{ROOT}/classified_pairs.tsv", f"{ROOT}/gate.json", f"{ROOT}/output_sha256.txt",
    resources:
        mem_mb=4000,
    shell:
        "{PY:q} {CODE:q} report --config {CONFIG:q}"


rule verification:
    input:
        rules.report_coi_origin_verification.output,
    default_target: True
