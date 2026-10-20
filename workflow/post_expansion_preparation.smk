"""Consume completed assemblies and prepare review inputs; mapping is unavailable."""

import csv
from pathlib import Path

configfile: "config/assembly_completion.json"

PREP = config["preparation_root"]
PY = config["python"]
with open(config["cohort"]) as handle:
    COHORT = list(csv.DictReader(handle, delimiter="\t"))
SOURCE_ROOTS = sorted({r["result_root"] for r in COHORT})
PRODUCT_SUFFIXES = {"assemblies": [".fasta"], "markers": [".gff"], "viruses": [".fna", ".tsv"], "checkv": [".tsv"], "bins": [".tar.gz", ".depth.tsv"], "checkm2": [".tsv", ".status"]}


rule stage_cohort_for_qc:
    input:
        configuration="config/assembly_completion.json",
        hold="config/post_expansion_remapping.json",
        samples=config["cohort"],
        code=["scripts/prepare_mapping_review.py", "scripts/assembly_completion.py"],
        completion=[f"{s}/stages/phase3_assembly.done" for s in SOURCE_ROOTS],
        summaries=[f"{s}/phase3_cohort/assembly_summary.tsv" for s in SOURCE_ROOTS],
        products=[f"{r['result_root']}/phase3_cohort/{kind}/{r['sample_id']}{suffix}" for r in COHORT for kind, suffixes in PRODUCT_SUFFIXES.items() for suffix in suffixes],
        records=[f"{r['input_root']}/provenance/{'phase3_cohort/inputs' if r['host_route'] != 'none' else 'trim'}/{r['sample_id']}.json" for r in COHORT],
    output:
        products=directory(f"{PREP}/phase3_cohort"),
        sources=f"{PREP}/assembly_sources.tsv",
        mapping=f"{PREP}/mapping_inputs.tsv",
        provenance=f"{PREP}/assembly_preparation.json",
    resources:
        mem_mb=2000,
        runtime=120,
    shell:
        "{PY:q} scripts/prepare_mapping_review.py stage --config {input.configuration:q}"


rule select_mags_for_qc:
    input:
        staged=rules.stage_cohort_for_qc.output,
        script="workflow/scripts/select_mags.py",
        config="config/phase3_catalog.json",
    output:
        genomes=directory(f"{PREP}/qc/candidates"),
        table=f"{PREP}/qc/candidates.tsv",
    shell:
        "{PY:q} {input.script:q} --membership {PREP:q}/phase3_cohort/membership.tsv --results {PREP:q}/phase3_cohort --config {input.config:q} --genomes {output.genomes:q} --table {output.table:q}"


rule report_catalog_qc_inputs:
    input:
        candidates=rules.select_mags_for_qc.output,
        sources=rules.stage_cohort_for_qc.output,
        hold="config/post_expansion_remapping.json",
        exclusions="config/catalog_domain_exclusions.tsv",
        original_audit="config/mag_domain_audit.json",
        old_catalog="data/results/phase3_catalog/mags/mag_catalog.tsv",
        code=["scripts/prepare_mapping_review.py", "scripts/assembly_completion.py"],
    output:
        candidates="figures/catalog_qc/candidate_review.tsv",
        contigs=f"{PREP}/qc/contig_review.tsv",
        report="figures/catalog_qc/preparation_summary.json",
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "{PY:q} scripts/prepare_mapping_review.py review"


rule qc_preparation:
    input: rules.report_catalog_qc_inputs.output
    default_target: True
