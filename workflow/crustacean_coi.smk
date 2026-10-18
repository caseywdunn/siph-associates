"""Isolated COI pilot using all retained host-depleted N. septata pairs."""

import csv
import json
from pathlib import Path

configfile: "config/crustacean_coi.json"

PANEL = config["panel"]
RUN = config["run"]
PY = config["python"]
DIAMOND = config["diamond"]
BLAST = config["blast_bin"]
BB = config["bbtools_bin"]
SEQKIT = config["seqkit"]
with open(config["manifest"]) as handle:
    SAMPLES = {r["analysis_sample_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
SOURCES = json.loads(Path(config["sources"]).read_text())


def reads(wildcards, mate):
    if wildcards.sample == "controls":
        return f"{RUN}/controls/R{mate}.fastq"
    return SAMPLES[wildcards.sample][f"r{mate}_path"]


def control_gate(wildcards):
    return [] if wildcards.sample == "controls" else [f"{RUN}/controls/PASSED"]


def expected_pairs(wildcards):
    count = 12 if wildcards.sample == "controls" else int(SAMPLES[wildcards.sample]["pairs"])
    return f"--expected-pairs {count}"

rule references:
    input:
        f"{PANEL}/target.masked.fna",
        f"{PANEL}/target.dmnd",
        f"{PANEL}/competitive.nsq",
        f"{PANEL}/coverage.tsv",

rule fetch_coi_references:
    input:
        sources=config["sources"],
        code="scripts/crustacean_coi/fetch_references.py",
    output:
        receipt=f"{config['reference_dir']}/downloads.json",
        files=[f"{config['reference_dir']}/{s['filename']}" for s in SOURCES],
    params:
        directory=config["reference_dir"],
    shell:
        "{PY:q} {input.code:q} --sources {input.sources:q} --output {params.directory:q}"

rule extract_coi:
    input:
        receipt=rules.fetch_coi_references.output.receipt,
        files=rules.fetch_coi_references.output.files,
        sources=config["sources"],
        code="scripts/crustacean_coi/prepare_reference.py",
    output:
        nt=f"{PANEL}/annotated.fna",
        aa=f"{PANEL}/annotated.faa",
        barcodes=f"{PANEL}/barcodes.fna",
        host=f"{PANEL}/host_barcodes.fna",
        annotated_metadata=f"{PANEL}/annotated_metadata.tsv",
        barcode_metadata=f"{PANEL}/barcode_metadata.tsv",
        summary=f"{PANEL}/extraction.json",
    params:
        directory=config["reference_dir"],
    resources:
        mem_mb=3000,
    benchmark:
        f"{PANEL}/benchmarks/extract.tsv"
    shell:
        "{PY:q} {input.code:q} extract --reference-dir {params.directory:q} "
        "--sources {input.sources:q} --output {PANEL:q}"

rule coi_frame_anchors:
    input:
        aa=f"{PANEL}/annotated.faa",
        metadata=f"{PANEL}/annotated_metadata.tsv",
        code="scripts/crustacean_coi/frame_anchors.py",
    output:
        f"{PANEL}/frame_anchors.faa",
    shell:
        "{PY:q} {input.code:q} {PANEL:q}"


rule annotated_coi_proteins:
    input:
        f"{PANEL}/frame_anchors.faa",
    output:
        f"{PANEL}/frame_anchors.dmnd",
    threads: 4
    resources:
        mem_mb=2000,
    shell:
        "{DIAMOND:q} makedb --in {input:q} --db {output:q} --threads {threads}"

rule validate_coi_barcode_frames:
    input:
        db=f"{PANEL}/frame_anchors.dmnd",
        query=f"{PANEL}/barcodes.fna",
    output:
        f"{PANEL}/barcode_frames.tsv",
    threads: 4
    resources:
        mem_mb=4000,
    benchmark:
        f"{PANEL}/benchmarks/barcode_frames.tsv"
    shell:
        "{DIAMOND:q} blastx --db {input.db:q} --query {input.query:q} "
        "--out {output:q} --query-gencode 5 --sensitive --min-orf 20 "
        "--max-target-seqs 1 --query-cover 80 --evalue 1e-5 --block-size 0.2 "
        "--index-chunks 1 --threads {threads} "
        "--outfmt 6 qseqid sseqid pident length qcovhsp qframe bitscore"

rule finalize_coi_panel:
    input:
        frames=f"{PANEL}/barcode_frames.tsv",
        extracted=rules.extract_coi.output,
        code="scripts/crustacean_coi/prepare_reference.py",
    output:
        nt=f"{PANEL}/target.fna",
        aa=f"{PANEL}/target.faa",
        competitive=f"{PANEL}/competitive.fna",
        metadata=f"{PANEL}/reference_metadata.tsv",
        coverage=f"{PANEL}/coverage.tsv",
        summary=f"{PANEL}/reference_summary.json",
    resources:
        mem_mb=3000,
    shell:
        "{PY:q} {input.code:q} finalize --output {PANEL:q}"



rule pilot:
    input:
        f"{RUN}/controls/PASSED",
        expand(f"{RUN}/{{sample}}/summary.json", sample=SAMPLES),


rule controls:
    input:
        f"{RUN}/controls/PASSED",


rule make_coi_controls:
    input:
        nt=f"{PANEL}/competitive.fna",
        metadata=f"{PANEL}/reference_metadata.tsv",
        code="scripts/crustacean_coi/controls.py",
    output:
        f"{RUN}/controls/R1.fastq",
        f"{RUN}/controls/R2.fastq",
        f"{RUN}/controls/truth.json",
    resources:
        mem_mb=2000,
    shell:
        "{PY:q} {input.code:q} generate --panel {PANEL:q} --output {RUN:q}/controls"


include: "rules/crustacean_coi_search.smk"


rule check_coi_controls:
    input:
        truth=f"{RUN}/controls/truth.json",
        summary=f"{RUN}/controls/summary.json",
        classified=f"{RUN}/controls/classified_pairs.tsv",
        code="scripts/crustacean_coi/controls.py",
    output:
        f"{RUN}/controls/control_results.json",
        f"{RUN}/controls/PASSED",
    shell:
        "{PY:q} {input.code:q} check --output {RUN:q}/controls"
