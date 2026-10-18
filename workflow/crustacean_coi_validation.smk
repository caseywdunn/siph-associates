"""Audit real COI candidates and test sensitivity with genera withheld.

Consume preserved pilot artifacts, build a derivative panel, reuse the pilot
search rules, and report recruitment separately from competitive assignment.
"""

import json
from pathlib import Path

configfile: "config/crustacean_coi_validation.json"

wildcard_constraints:
    sample="[^/]+",

RUN = config["validation_run"]
PANEL = f"{RUN}/panel"
PY, DIAMOND = config["python"], config["diamond"]
BLAST, BB, SEQKIT = config["blast_bin"], config["bbtools_bin"], config["seqkit"]
SOURCE = config["source_panel"]
PILOT = config["pilot_run"]
CONFIG = "config/crustacean_coi_validation.json"
config["manifest"] = f"{RUN}/heldout/truth.json"
# Track the known sample list without requiring generated inputs at parse time.
import csv
with open("config/crustacean_coi_pilot.tsv") as handle:
    SAMPLES = [r["analysis_sample_id"] for r in csv.DictReader(handle, delimiter="\t")]


def reads(wildcards, mate):
    return f"{RUN}/heldout/R{mate}.fastq"


def control_gate(wildcards):
    return []


def expected_pairs(wildcards):
    return f"--expected-pairs-json {RUN}/heldout/truth.json"


rule prepare_heldout_coi_panel:
    input:
        [f"{SOURCE}/{name}" for name in ("target.fna", "target.faa", "competitive.fna", "reference_metadata.tsv")],
        [f"{PILOT}/controls/{name}" for name in ("R1.fastq", "R2.fastq", "truth.json", "PASSED")],
        "scripts/crustacean_coi/validation.py", CONFIG,
    output:
        [f"{PANEL}/{name}" for name in ("target.fna", "target.faa", "competitive.fna", "reference_metadata.tsv")],
        [f"{RUN}/heldout/{name}" for name in ("R1.fastq", "R2.fastq", "truth.json")],
        f"{RUN}/design.json",
    resources:
        mem_mb=2000,
    shell:
        "{PY:q} scripts/crustacean_coi/validation.py prepare --config {CONFIG:q}"


include: "rules/crustacean_coi_search.smk"


rule prepare_pilot_audit:
    input:
        expand(f"{PILOT}/{{sample}}/candidates.fna", sample=SAMPLES),
        expand(f"{PILOT}/{{sample}}/candidate_pairs.tsv", sample=SAMPLES),
        f"{SOURCE}/target.fna", "scripts/crustacean_coi/validation.py", CONFIG,
    output:
        f"{RUN}/audit/queries.fna", f"{RUN}/audit/query_pairs.tsv",
    shell:
        "{PY:q} scripts/crustacean_coi/validation.py audit --config {CONFIG:q}"


rule audit_whole_mitochondria_blastn:
    input:
        query=f"{RUN}/audit/queries.fna",
        db=config["mito_db"] + ".nsq",
    output:
        f"{RUN}/audit/mitochondrial.tsv",
    log:
        f"{RUN}/audit/mitochondrial.log",
    threads: 4
    resources:
        mem_mb=3000,
    shell:
        "{BLAST:q}/blastn -query {input.query:q} -db {config[mito_db]:q} "
        "-task blastn -word_size 7 -dust yes -evalue 1e-5 -max_hsps 1 "
        "-max_target_seqs 100 -num_threads {threads} -out {output:q} "
        "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle' 2> {log:q}"


rule audit_swissprot_diamond:
    input:
        query=f"{RUN}/audit/queries.fna", db=config["swissprot_db"],
    output:
        f"{RUN}/audit/swissprot_code{{code}}.tsv",
    wildcard_constraints:
        code="[15]",
    log:
        f"{RUN}/audit/swissprot_code{{code}}.log",
    threads: 4
    resources:
        mem_mb=4000,
    shell:
        "{DIAMOND:q} blastx --query {input.query:q} --db {input.db:q} --out {output:q} "
        "--query-gencode {wildcards.code} --very-sensitive --min-orf 20 --evalue 1e-3 "
        "--max-target-seqs 100 --block-size 0.2 --index-chunks 1 --threads {threads} "
        "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle "
        "2> {log:q}"


rule revise_pilot_classification:
    input:
        alignment=f"{PILOT}/{{sample}}/competitive.tsv",
        pairs=f"{PILOT}/{{sample}}/candidate_pairs.tsv",
        counts=f"{PILOT}/{{sample}}/recruitment.json",
        metadata=f"{SOURCE}/reference_metadata.tsv",
        code="scripts/crustacean_coi/summarize.py",
    output:
        f"{RUN}/reclassified/{{sample}}/summary.json",
        f"{RUN}/reclassified/{{sample}}/classified_pairs.tsv",
    params:
        source=lambda w: f"{PILOT}/{w.sample}",
        destination=lambda w: f"{RUN}/reclassified/{w.sample}",
    shell:
        "{PY:q} {input.code:q} --directory {params.source:q} --metadata {input.metadata:q} --output {params.destination:q}"


rule report_coi_validation:
    input:
        f"{RUN}/heldout/summary.json", f"{RUN}/heldout/classified_pairs.tsv",
        f"{RUN}/heldout/truth.json", f"{RUN}/audit/mitochondrial.tsv",
        expand(f"{RUN}/audit/swissprot_code{{code}}.tsv", code=[1, 5]),
        expand(f"{RUN}/reclassified/{{sample}}/summary.json", sample=SAMPLES + ["controls"]),
        "scripts/crustacean_coi/validation.py", CONFIG,
    output:
        f"{RUN}/sensitivity_summary.json", f"{RUN}/sensitivity_pairs.tsv",
        f"{RUN}/output_sha256.txt",
    shell:
        "{PY:q} scripts/crustacean_coi/validation.py report --config {CONFIG:q}"


include: "rules/crustacean_coi_protein_validation.smk"



rule audit_bacterial_catalog_bwa:
    input:
        query=f"{RUN}/audit/queries.fna",
        genome="data/results/phase4_catalog/v2/bacterial_catalog.fna",
        index="data/results/phase4_catalog/v2/bacterial_catalog.fna.bwt",
        metadata="data/results/phase4_catalog/v2/bacterial_catalog.manifest.tsv",
    output:
        f"{RUN}/audit/bacterial_catalog.sam",
    log:
        f"{RUN}/audit/bacterial_catalog.log",
    threads: 4
    resources:
        mem_mb=7000,
    shell:
        "/gpfs/gibbs/project/dunn/cwd7/conda_envs/assembly_env/bin/bwa mem "
        "-a -t {threads} {input.genome:q} {input.query:q} > {output:q} 2> {log:q}"

rule audit_na35_assembly_blastn:
    input:
        query=f"{RUN}/audit/queries.fna",
        assembly="data/results/phase3_cohort/assemblies/Ahuja2026__NA35.fasta",
    output:
        f"{RUN}/audit/na35_assembly.tsv",
    shell:
        "{BLAST:q}/blastn -query {input.query:q} -subject {input.assembly:q} "
        "-task blastn -word_size 11 -dust yes -evalue 1e-8 -max_hsps 1 "
        "-out {output:q} -outfmt '6 qseqid sseqid pident length qstart qend sstart send bitscore evalue'"

rule report_local_origin_audit:
    input:
        rules.audit_bacterial_catalog_bwa.output,
        f"{RUN}/reclassified/Ahuja2026__NA35/classified_pairs.tsv",
        "data/results/phase4_catalog/v2/bacterial_catalog.manifest.tsv",
        "scripts/crustacean_coi/summarize_local_audit.py",
    output:
        f"{RUN}/audit/na35_bacterial_matches.tsv",
        f"{RUN}/audit/na35_audited_pairs.tsv",
    shell:
        "{PY:q} scripts/crustacean_coi/summarize_local_audit.py --root {RUN:q} "
        "--catalog data/results/phase4_catalog/v2/bacterial_catalog.manifest.tsv"


rule local_origin_validation:
    input:
        rules.audit_bacterial_catalog_bwa.output,
        rules.audit_na35_assembly_blastn.output,
        rules.report_local_origin_audit.output,


rule validation:
    input:
        rules.report_coi_validation.output,
    default_target: True
