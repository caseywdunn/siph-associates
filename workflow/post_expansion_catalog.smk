"""Curate viruses and compare genomes before freezing mapping targets.

Reuse completed per-library products and shared catalog methods. No read mapping
rules are included. MAG inclusion follows the separate domain review.
"""

import json
from pathlib import Path

configfile: "config/post_expansion_catalog.json"

WORK = config["work_root"]
SCRATCH = config["scratch_root"]
PHASE3_CATALOG_CONFIG_PATH = Path("config/phase3_catalog.json")
PHASE3_CATALOG = json.loads(PHASE3_CATALOG_CONFIG_PATH.read_text())
PHASE3_ASSEMBLY_MEMBERSHIP_PATH = Path(WORK) / "phase3_cohort/membership.tsv"
PHASE4_CATALOG = {"catalog_version": config["catalog_version"]}
for resources in PHASE3_CATALOG["resources"].values():
    resources["partition"] = "ycga"

def catalog_resources(name):
    return {k: v for k, v in PHASE3_CATALOG["resources"][name].items() if k != "threads"}

include: "rules/phase3_catalog.smk"

rule classify_new_candidates_gtdbtk:
    input:
        batch="config/catalog_qc_new_genomes.tsv",
        fasta=[line.split("\t")[0] for line in Path("config/catalog_qc_new_genomes.tsv").read_text().splitlines()],
        reference=ancient(f"{SCRATCH}/phase3_catalog/gtdbtk_r220"),
    output:
        directory(f"{WORK}/qc/new_gtdb"),
    params:
        tools=CATALOG_TOOLS["gtdbtk"],
        mash=f"{SCRATCH}/phase3_catalog/new_gtdb.msh",
    log: f"{WORK}/logs/qc/new_gtdb.log"
    threads: 16
    resources: mem_mb=150000, runtime=720, partition="ycga"
    shell:
        r"""
        export PATH={params.tools:q}:$PATH GTDBTK_DATA_PATH={input.reference:q}
        gtdbtk classify_wf --batchfile {input.batch:q} --out_dir {output:q} \
          --cpus {threads} --pplacer_cpus 8 --mash_db {params.mash:q} > {log:q} 2>&1
        """

rule compare_review_candidates_skani:
    input:
        genomes=f"{WORK}/qc/candidates",
        table=f"{WORK}/qc/candidates.tsv",
    output: f"{WORK}/qc/candidates.skani.tsv"
    params: tools=CATALOG_TOOLS["skani"]
    log: f"{WORK}/logs/qc/candidates_skani.log"
    threads: 8
    resources: mem_mb=16000, runtime=60, partition="ycga"
    shell:
        "{params.tools:q}/skani triangle -t {threads} -E {input.genomes:q}/*.fa > {output:q}.tmp 2> {log:q}; mv {output:q}.tmp {output:q}"

rule catalog_evidence:
    input:
        expand(f"{WORK}/phase3_catalog/votus/{{virus_set}}.votus.tsv", virus_set=["associate", "endogenous_candidate"]),
        rules.classify_new_candidates_gtdbtk.output,
        rules.compare_review_candidates_skani.output,
    default_target: True
