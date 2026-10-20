"""Retain gene-level viral evidence and resolve cellular flanks before mapping."""

import json
from pathlib import Path

configfile: "config/post_expansion_catalog.json"
WORK = config["work_root"]
QC = f"{WORK}/qc/viruses"
ASSEMBLY = json.loads(Path("config/phase3_assembly.json").read_text())
TOOLS = ASSEMBLY["tool_prefixes"]["genomad"] + "/bin"

rule annotate_viral_candidates_genomad:
    input:
        fasta=f"{WORK}/phase3_catalog/viruses/included.fna",
        database=ASSEMBLY["databases"]["genomad"],
    output: directory(f"{QC}/genomad")
    threads: 16
    resources: mem_mb=32000, runtime=360, partition="ycga"
    log: f"{WORK}/logs/qc/viral_genomad.log"
    shell:
        "export PATH={TOOLS:q}:$PATH; genomad annotate --threads {threads} {input.fasta:q} {output:q} {input.database:q} > {log:q} 2>&1"

rule assess_viral_boundaries_checkv:
    input:
        fasta=f"{WORK}/phase3_catalog/viruses/included.fna",
        database=ASSEMBLY["databases"]["checkv"],
    output: directory(f"{QC}/checkv")
    threads: 16
    resources: mem_mb=32000, runtime=240, partition="ycga"
    log: f"{WORK}/logs/qc/viral_checkv.log"
    shell:
        "export PATH={TOOLS:q}:$PATH; checkv end_to_end {input.fasta:q} {output:q} -d {input.database:q} -t {threads} > {log:q} 2>&1"

rule viral_evidence:
    input:
        rules.annotate_viral_candidates_genomad.output,
        rules.assess_viral_boundaries_checkv.output,
    default_target: True
