"""Adjudicate discordant domain evidence with nr; no read mapping is included."""

configfile: "config/catalog_nr_review.json"
ROOT = config["run"]
PYTHON = config["python"]
DIAMOND = config["diamond"]
CONFIG = config["settings_path"]
SCRIPT = "scripts/mag_domain_audit"
REF = config["reference_root"]
OUTPUT_FIELDS = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen stitle"

module domain:
    snakefile: "mag_domain_audit.smk"
    config: config

rule select_discordant_genes:
    input:
        bins="figures/catalog_qc/domain_review.tsv",
        contigs="data/results/post_expansion/2026-10-20/qc/contig_evidence.tsv",
        genes=config["screen_root"] + "/genes.faa",
        metadata=config["screen_root"] + "/genes.tsv",
        assignments=config["screen_root"] + "/review/genes_swissprot.tsv",
        markers=config["screen_root"] + "/archaeal_marker_hits.tsv",
        code=["scripts/select_catalog_nr_queries.py", "scripts/assembly_completion.py"],
        configuration=CONFIG,
    output:
        fasta=f"{ROOT}/genes.faa",
        metadata=f"{ROOT}/genes.tsv",
        markers=f"{ROOT}/archaeal_marker_hits.tsv",
        selection=f"{ROOT}/query_selection.tsv",
        provenance=f"{ROOT}/selection.json",
    resources: mem_mb=4000, runtime=60, partition="ycga"
    shell: "{PYTHON:q} scripts/select_catalog_nr_queries.py"

use rule search_mag_genes_nr_diamond from domain as search_mag_genes_nr_diamond with:
    input:
        query=f"{ROOT}/genes.faa",
        database=config["nr_database"],
        reference_hash=f"{REF}/nr_sha256.txt",
        controls=config["screen_root"] + "/controls/CONTROLS_PASSED",
    threads: 32
    resources: mem_mb=96000, runtime=1440, partition="ycga"

use rule report_audit_assignments from domain as report_audit_assignments with:
    resources: mem_mb=12000, runtime=120, partition="ycga"

rule nr_review:
    input: f"{ROOT}/review/genes_nr.counts.tsv"
    default_target: True
