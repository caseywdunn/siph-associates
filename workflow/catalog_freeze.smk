"""Curate, dereplicate and index the expanded catalogs; never map library reads."""

import json
from pathlib import Path

configfile: "config/catalog_freeze.json"
WORK, SCRATCH = config["work_root"], config["scratch_root"]
PREP = config["preparation_root"]
PY = "/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python"
PHASE3_CATALOG_CONFIG_PATH = Path("config/phase3_catalog.json")
PHASE4_CATALOG_CONFIG_PATH = Path("config/catalog_freeze_policy.json")
PHASE3_CATALOG = json.loads(PHASE3_CATALOG_CONFIG_PATH.read_text())
PHASE4_CATALOG = json.loads(PHASE4_CATALOG_CONFIG_PATH.read_text())
# Snakemake 7 formats imported shell expressions in the outer namespace.
CATALOG = f"{WORK}/phase3_catalog"
CATALOG_SCRATCH = f"{SCRATCH}/phase3_catalog"
CATALOG_TOOLS = {k: v + "/bin" for k, v in PHASE3_CATALOG["tool_prefixes"].items()}
CATALOG_RESOURCES = PHASE3_CATALOG["resources"]
VIRUS_SETS = {"associate": "VOTU", "endogenous_candidate": "EVEC"}
HOST_ROUTES = sorted(config["databases"]["host_references"])
P4 = f"{WORK}/phase4_catalog"
P4_VERSION = PHASE4_CATALOG["catalog_version"]
P4_SCRATCH = f"{SCRATCH}/phase4_catalog"
P4_BACTERIAL = f"{P4}/{P4_VERSION}/bacterial_catalog"
P4_TOOLS = {k: v + "/bin" for k, v in PHASE4_CATALOG["tool_prefixes"].items()}
P4_VIRAL_SETS = ("associate", "endogenous_candidate")

module shared:
    snakefile: "catalog_shared.smk"
    config: config

rule merge_reference_evidence:
    input:
        code=["scripts/prepare_expanded_references.py", "scripts/assembly_completion.py"],
        old="data/results/phase4_catalog/references",
        new=f"{PREP}/phase4_catalog/references/new",
        gtdb=f"{PREP}/phase4_catalog/references/new/gtdb",
        checkm=f"{PREP}/phase4_catalog/references/new/checkm2_quality.tsv",
        nominated=f"{PREP}/phase4_catalog/references/nominated.tsv",
        decoys=f"{PREP}/phase4_catalog/references/decoy_candidates.tsv",
    output:
        tables=expand(f"{P4}/references/{{name}}", name=["nominated.tsv", "decoy_candidates.tsv", "staged.tsv", "ncbi_status.tsv", "substitutes.tsv", "checkm2_quality.tsv", "oceandna_gtdbtk.bac120.summary.tsv", "oceandna_gtdbtk.ar53.summary.tsv", "qc_sources.tsv"]),
        genomes=directory(f"{P4}/references/genomes"),
        substitutes=directory(f"{P4}/references/substitutes"),
    params: destination=f"{P4}/references"
    resources: mem_mb=4000, runtime=120, partition="ycga"
    shell: "{PY:q} scripts/prepare_expanded_references.py merge --destination {params.destination:q}"

rule prepare_reviewed_mags:
    input:
        code=["scripts/curate_catalog_inputs.py", "scripts/assembly_completion.py"],
        decisions="figures/catalog_qc/bin_decisions.tsv",
        contigs="figures/catalog_qc/contig_decisions.tsv",
        candidates=f"{PREP}/qc/candidates.tsv",
        genomes=f"{PREP}/qc/candidates",
        previous="figures/catalog_qc/candidate_review.tsv",
        old_taxonomy=expand("data/results/phase3_catalog/mags/gtdbtk.{domain}.summary.tsv", domain=["bac120", "ar53"]),
        new_taxonomy=f"{PREP}/qc/new_gtdb/gtdbtk.bac120.summary.tsv",
    output:
        genomes=directory(f"{CATALOG}/mags/candidates"),
        changed=directory(f"{CATALOG}/mags/changed"),
        table=f"{CATALOG}/mags/before_quality.tsv",
        taxonomy=expand(f"{CATALOG}/mags/gtdbtk.{{domain}}.summary.tsv", domain=["bac120", "ar53"]),
    resources: mem_mb=4000, runtime=60, partition="ycga"
    shell: "{PY:q} scripts/curate_catalog_inputs.py prepare"

rule reassess_curated_mags_checkm2:
    input:
        genomes=rules.prepare_reviewed_mags.output.changed,
        database=PHASE4_CATALOG["databases"]["checkm2"],
    output: directory(f"{CATALOG}/mags/changed_checkm2")
    params: tools=P4_TOOLS["checkm2"]
    threads: 16
    resources: mem_mb=64000, runtime=480, partition="ycga"
    log: f"{WORK}/logs/curated_checkm2.log"
    shell:
        "export PATH={params.tools:q}:$PATH; checkm2 predict --threads {threads} --input {input.genomes:q} --extension fa --database_path {input.database:q} --output-directory {output:q} --force > {log:q} 2>&1"

rule accept_curated_mag_quality:
    input:
        previous=rules.prepare_reviewed_mags.output.table,
        checkm2=rules.reassess_curated_mags_checkm2.output,
        code=["scripts/curate_catalog_inputs.py", "scripts/assembly_completion.py"],
        policy="config/phase3_catalog.json",
    output:
        table=f"{CATALOG}/mags/candidates.tsv",
        report="figures/catalog_qc/curated_mag_quality.tsv",
    resources: mem_mb=4000, runtime=60, partition="ycga"
    shell: "{PY:q} scripts/curate_catalog_inputs.py quality"

rule accept_curated_viral_regions:
    input:
        regions=f"{PREP}/qc/curation/viral_regions.fna",
        table=f"{PREP}/qc/curation/viral_regions.tsv",
        quality=f"{PREP}/qc/curation/viral_checkv/quality_summary.tsv",
        boundaries=f"{PREP}/qc/curation/viral_checkv/contamination.tsv",
        code=["scripts/curate_catalog_inputs.py", "scripts/assembly_completion.py"],
        policy="config/phase3_catalog.json",
    output:
        table=f"{CATALOG}/viruses/included.tsv",
        fasta=f"{CATALOG}/viruses/included.fna",
        review="figures/catalog_qc/viral_region_decisions.tsv",
    resources: mem_mb=4000, runtime=60, partition="ycga"
    shell: "{PY:q} scripts/curate_catalog_inputs.py viruses"

use rule compare_mags_skani, build_mag_catalog, align_viruses_host_minimap2, flag_endogenous_viruses, cluster_votus_checkv, build_votu_catalog from shared as *
use rule compare_catalog_genomes_skani, build_bacterial_catalog, index_bacterial_catalog_bwa, build_viral_catalog, call_crispr_spacers_minced, match_spacers_votus_blastn, validate_phase4_catalog from shared as *

rule index_combined_viral_catalog:
    input:
        expand(f"{P4}/{P4_VERSION}/viral_{{kind}}.fna", kind=P4_VIRAL_SETS),
    output:
        fasta=f"{P4}/{P4_VERSION}/viral_catalog.fna",
        index=multiext(f"{P4}/{P4_VERSION}/viral_catalog.fna", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"),
    params: tools=P4_TOOLS["bwa"]
    threads: 1
    resources: mem_mb=8000, runtime=120, partition="ycga"
    log: f"{WORK}/logs/viral_index.log"
    shell:
        "export PATH={params.tools:q}:$PATH; cat {input:q} > {output.fasta:q}; bwa index {output.fasta:q} > {log:q} 2>&1; samtools faidx {output.fasta:q} >> {log:q} 2>&1"

rule record_catalog_readiness:
    input:
        validation=rules.validate_phase4_catalog.output,
        viral=rules.index_combined_viral_catalog.output,
        bacterial=rules.index_bacterial_catalog_bwa.output,
        mags=rules.build_mag_catalog.output.catalog,
        decisions=["figures/catalog_qc/bin_decisions.tsv", "figures/catalog_qc/contig_decisions.tsv", "figures/catalog_qc/curated_mag_quality.tsv", "figures/catalog_qc/viral_region_decisions.tsv"],
        code=["scripts/record_catalog_ready.py", "scripts/assembly_completion.py"],
        configuration="config/catalog_freeze.json",
        hold="config/post_expansion_remapping.json",
    output:
        report="figures/catalog_qc/catalog_v3/catalog_ready.json",
        inventory="figures/catalog_qc/catalog_v3/inventory.tsv",
        crosswalk="figures/catalog_qc/catalog_v3/catalog_id_crosswalk.tsv",
    resources: mem_mb=4000, runtime=60, partition="ycga"
    shell: "{PY:q} scripts/record_catalog_ready.py"

rule catalog_ready:
    input: rules.record_catalog_readiness.output
    default_target: True
