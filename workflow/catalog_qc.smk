"""Audit every candidate using the shared gene/translated-DNA domain searches.

Consume completed candidates and frozen reference databases. This bounded entry
point imports no reference downloads, assemblies, or cohort read-mapping rules.
"""

configfile: "config/catalog_qc.json"

# Snakemake 7 formats imported shell expressions in the importing namespace.
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

use rule prepare_audit_inputs, collect_audit_genes, search_audit_controls_diamond, validate_audit_controls from domain as *

use rule identify_mag_markers_gtdbtk from domain as identify_mag_markers_gtdbtk with:
    threads: 16
    resources: mem_mb=32000, runtime=240, partition="ycga"

use rule search_mag_genes_swissprot_diamond from domain as search_mag_genes_swissprot_diamond with:
    threads: 16
    resources: mem_mb=16000, runtime=240, partition="ycga"

use rule search_mag_tiles_swissprot_diamond from domain as search_mag_tiles_swissprot_diamond with:
    threads: 16
    resources: mem_mb=16000, runtime=240, partition="ycga"

use rule report_audit_assignments from domain as report_audit_assignments with:
    resources: mem_mb=48000, runtime=120, partition="ycga"

rule summarize_domain_evidence:
    input:
        reports=expand(config["run"] + "/review/{review}.tsv", review=["genes_swissprot", "tiles_swissprot"]),
        markers=config["run"] + "/gtdb_identify",
        candidates="figures/catalog_qc/candidate_review.tsv",
        contigs="data/results/post_expansion/2026-10-20/qc/contig_review.tsv",
        config="config/catalog_qc.json",
        code=["scripts/summarize_catalog_qc.py", "scripts/assembly_completion.py"],
    output:
        bins="figures/catalog_qc/domain_review.tsv",
        contigs="data/results/post_expansion/2026-10-20/qc/contig_evidence.tsv",
    resources: mem_mb=4000, runtime=60, partition="ycga"
    shell:
        "{PYTHON:q} scripts/summarize_catalog_qc.py"

rule domain_screen:
    input:
        expand(config["run"] + "/review/{review}.counts.tsv", review=["genes_swissprot", "tiles_swissprot"]),
        config["run"] + "/controls/CONTROLS_PASSED",
        rules.summarize_domain_evidence.output,
    default_target: True
