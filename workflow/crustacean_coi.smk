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

rule mask_coi_seeds:
    input:
        f"{PANEL}/target.fna",
    output:
        f"{PANEL}/target.masked.fna",
    resources:
        mem_mb=2000,
    shell:
        "{BLAST:q}/dustmasker -in {input:q} -outfmt fasta | "
        "{SEQKIT:q} replace --threads 1 -s -p '[a-z]' -r N > {output:q}"

rule index_target_coi_proteins:
    input:
        f"{PANEL}/target.faa",
    output:
        f"{PANEL}/target.dmnd",
    threads: 4
    resources:
        mem_mb=2000,
    shell:
        "{DIAMOND:q} makedb --in {input:q} --db {output:q} --threads {threads}"

rule index_competitive_coi:
    input:
        f"{PANEL}/competitive.fna",
    output:
        f"{PANEL}/competitive.nsq",
        f"{PANEL}/competitive.nhr",
        f"{PANEL}/competitive.nin",
    shell:
        "{BLAST:q}/makeblastdb -in {input:q} -dbtype nucl -parse_seqids "
        "-blastdb_version 4 -out {PANEL:q}/competitive"


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


rule nucleotide_coi_recruitment:
    input:
        r1=lambda w: reads(w, 1),
        r2=lambda w: reads(w, 2),
        reference=f"{PANEL}/target.masked.fna",
        controls=control_gate,
        manifest=config["manifest"],
    output:
        r1=f"{RUN}/{{sample}}/bbduk_R1.fastq",
        r2=f"{RUN}/{{sample}}/bbduk_R2.fastq",
        stats=f"{RUN}/{{sample}}/bbduk_stats.tsv",
    log:
        f"{RUN}/{{sample}}/bbduk.log",
    threads: 4
    resources:
        mem_mb=9000,
    benchmark:
        f"{RUN}/{{sample}}/benchmarks/bbduk.tsv"
    shell:
        "{BB:q}/bbduk.sh -Xmx8g in={input.r1:q} in2={input.r2:q} "
        "ref={input.reference:q} outm={output.r1:q} outm2={output.r2:q} "
        "stats={output.stats:q} k=25 hdist=0 maskmiddle=f rcomp=t "
        "rieb=t qtrim=f minlen=0 overwrite=t threads={threads} 2> {log:q}"


rule translated_coi_recruitment:
    input:
        reads=lambda w: reads(w, w.mate),
        db=f"{PANEL}/target.dmnd",
        controls=control_gate,
        manifest=config["manifest"],
    output:
        f"{RUN}/{{sample}}/diamond_R{{mate}}.tsv",
    wildcard_constraints:
        mate="[12]",
    log:
        f"{RUN}/{{sample}}/diamond_R{{mate}}.log",
    threads: 4
    resources:
        mem_mb=5000,
    benchmark:
        f"{RUN}/{{sample}}/benchmarks/diamond_R{{mate}}.tsv"
    shell:
        "{DIAMOND:q} blastx --query {input.reads:q} --db {input.db:q} --out {output:q} "
        "--query-gencode 5 --sensitive --min-orf 20 --evalue 1e-3 --query-cover 60 "
        "--max-target-seqs 1 --block-size 0.2 --index-chunks 1 --threads {threads} "
        "--outfmt 6 qseqid sseqid pident length bitscore evalue 2> {log:q}"


rule collect_coi_candidates:
    input:
        r1=lambda w: reads(w, 1),
        r2=lambda w: reads(w, 2),
        bbduk=f"{RUN}/{{sample}}/bbduk_R1.fastq",
        diamond1=f"{RUN}/{{sample}}/diamond_R1.tsv",
        diamond2=f"{RUN}/{{sample}}/diamond_R2.tsv",
        code="scripts/crustacean_coi/collect_candidates.py",
    output:
        fasta=f"{RUN}/{{sample}}/candidates.fna",
        pairs=f"{RUN}/{{sample}}/candidate_pairs.tsv",
        metrics=f"{RUN}/{{sample}}/recruitment.json",
    params:
        pairs=lambda w: 12 if w.sample == "controls" else int(SAMPLES[w.sample]["pairs"]),
        directory=lambda w: f"{RUN}/{w.sample}",
    resources:
        mem_mb=2000,
    benchmark:
        f"{RUN}/{{sample}}/benchmarks/collect.tsv"
    shell:
        "{PY:q} {input.code:q} --r1 {input.r1:q} --r2 {input.r2:q} "
        "--bbduk {input.bbduk:q} --diamond {input.diamond1:q} {input.diamond2:q} "
        "--expected-pairs {params.pairs} --output {params.directory:q}"


rule competitive_coi_alignment:
    input:
        query=f"{RUN}/{{sample}}/candidates.fna",
        db=rules.index_competitive_coi.output,
    output:
        f"{RUN}/{{sample}}/competitive.tsv",
    log:
        f"{RUN}/{{sample}}/competitive.log",
    threads: 4
    resources:
        mem_mb=2000,
    benchmark:
        f"{RUN}/{{sample}}/benchmarks/competitive.tsv"
    shell:
        "if test -s {input.query:q}; then "
        "{BLAST:q}/blastn -query {input.query:q} -db {PANEL:q}/competitive "
        "-task blastn -word_size 11 -dust yes -evalue 1e-8 -max_hsps 1 "
        "-max_target_seqs 200 -num_threads {threads} -out {output:q} "
        "-outfmt '6 qseqid sseqid pident length qstart qend sstart send bitscore evalue' "
        "2> {log:q}; else touch {output:q} {log:q}; fi"


rule summarize_coi:
    input:
        alignment=f"{RUN}/{{sample}}/competitive.tsv",
        pairs=f"{RUN}/{{sample}}/candidate_pairs.tsv",
        metrics=f"{RUN}/{{sample}}/recruitment.json",
        metadata=f"{PANEL}/reference_metadata.tsv",
        code="scripts/crustacean_coi/summarize.py",
    output:
        f"{RUN}/{{sample}}/classified_pairs.tsv",
        f"{RUN}/{{sample}}/summary.json",
    params:
        directory=lambda w: f"{RUN}/{w.sample}",
    shell:
        "{PY:q} {input.code:q} --directory {params.directory:q} --metadata {input.metadata:q}"


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
