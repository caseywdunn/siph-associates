"""Resolve MAG host fragments and viral boundaries before catalog finalization.

These are assembly-to-reference comparisons, not library read mapping.
"""

import json
from pathlib import Path

configfile: "config/post_expansion_catalog.json"
WORK = config["work_root"]
QC = f"{WORK}/qc/curation"
PY = "/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python"
TOOLS = config["tool_prefixes"]["assembly"] + "/bin"
ASSEMBLY = json.loads(Path("config/phase3_assembly.json").read_text())
REFS = dict(config["databases"]["host_references"], Perkinsela=f"{QC}/perkinsela_controls.fna")

rule prepare_curation_inputs:
    input:
        code=["scripts/prepare_catalog_curation.py", "scripts/assembly_completion.py"],
        bins="figures/catalog_qc/domain_review.tsv",
        genomes=f"{WORK}/qc/candidates",
        controls="config/mag_domain_audit.json",
        viruses=f"{WORK}/phase3_catalog/viruses/included.fna",
        table=f"{WORK}/phase3_catalog/viruses/included.tsv",
        boundaries=f"{WORK}/qc/viruses/checkv/contamination.tsv",
        genes=f"{WORK}/qc/viruses/genomad/included_annotate/included_genes.tsv",
    output:
        bacteria=f"{QC}/bacterial_candidates.fna",
        controls=f"{QC}/perkinsela_controls.fna",
        viruses=f"{QC}/viral_regions.fna",
        table=f"{QC}/viral_regions.tsv",
    resources: mem_mb=4000, runtime=60, partition="ycga"
    shell: "{PY:q} scripts/prepare_catalog_curation.py"

rule align_candidate_contigs_minimap2:
    input:
        contigs=rules.prepare_curation_inputs.output.bacteria,
        reference=lambda w: REFS[w.reference],
    output: f"{QC}/candidate_vs_{{reference}}.paf"
    threads: 8
    resources:
        mem_mb=lambda w: 64000 if w.reference == "H_rubrum" else 32000,
        runtime=120, partition="ycga"
    log: f"{WORK}/logs/curation/candidate_vs_{{reference}}.log"
    shell:
        "{TOOLS:q}/minimap2 -x asm20 -t {threads} --secondary=no -c {input.reference:q} {input.contigs:q} > {output:q}.tmp 2> {log:q}; mv {output:q}.tmp {output:q}"

rule reassess_trimmed_viruses_checkv:
    input:
        fasta=rules.prepare_curation_inputs.output.viruses,
        database=ASSEMBLY["databases"]["checkv"],
    output: directory(f"{QC}/viral_checkv")
    params: tools=ASSEMBLY["tool_prefixes"]["genomad"] + "/bin"
    threads: 16
    resources: mem_mb=32000, runtime=240, partition="ycga"
    log: f"{WORK}/logs/curation/viral_checkv.log"
    shell:
        "export PATH={params.tools:q}:$PATH; checkv end_to_end {input.fasta:q} {output:q} -d {input.database:q} -t {threads} > {log:q} 2>&1"

rule curation_evidence:
    input:
        expand(f"{QC}/candidate_vs_{{reference}}.paf", reference=sorted(REFS)),
        rules.reassess_trimmed_viruses_checkv.output,
    default_target: True
