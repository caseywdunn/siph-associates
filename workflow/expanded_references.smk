"""Nominate references from all four studies and QC newly nominated sequences.

Completed screen tables and reference QC are reused. There are no FASTQ
processing or cohort mapping rules in this entry point.
"""

import csv
import json
from pathlib import Path

configfile: "config/post_expansion_catalog.json"
WORK = config["work_root"]
REF = f"{WORK}/phase4_catalog/references"
NEW = f"{REF}/new"
BARG = "data/results/bargmannia/2026-10-19"
PY = "/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python"
P4 = json.loads(Path("config/phase4_catalog.json").read_text())
SAMPLES = [r["analysis_sample_id"] for r in csv.DictReader(open("config/bargmannia_samples.tsv"), delimiter="\t")]

rule extend_screen_nominations:
    input:
        original="data/results/aggregation/cohort/candidate_nominations.tsv",
        samples="config/bargmannia_samples.tsv",
        code=["scripts/extend_reference_nominations.py", "scripts/assembly_completion.py", "workflow/scripts/screen_nominations.py"],
        screens=[f"{BARG}/screens/{kind}/{s}{suffix}" for s in SAMPLES for kind, suffix in [("bracken", ".G.tsv"), ("bracken", ".S.tsv"), ("sylph", ".profile.tsv"), ("phyloflash", ".tar.gz")]],
        records=[f"{BARG}/provenance/{kind}/{s}.json" for s in SAMPLES for kind in ["trim", "kraken_bracken", "sylph", "phyloflash"]],
    output:
        table=f"{WORK}/aggregation/cohort/candidate_nominations.tsv",
        provenance=f"{WORK}/aggregation/cohort/provenance.json",
    resources: mem_mb=4000, runtime=120, partition="ycga"
    shell: "{PY:q} scripts/extend_reference_nominations.py"

rule nominate_expanded_references:
    input:
        nominations=rules.extend_screen_nominations.output.table,
        taxonomy=P4["references"]["gtdb_taxonomy"],
        config="config/phase4_catalog.json",
        code="workflow/scripts/nominate_references.py",
    output:
        references=f"{REF}/nominated.tsv",
        decoys=f"{REF}/decoy_candidates.redrawn.tsv",
    resources: mem_mb=4000, runtime=60, partition="ycga"
    shell:
        "{PY:q} {input.code:q} --nominations {input.nominations:q} --taxonomy {input.taxonomy:q} --config {input.config:q} --references {output.references:q} --decoys {output.decoys:q}"

rule select_new_references:
    input:
        nominations=rules.nominate_expanded_references.output.references,
        prior="data/results/phase4_catalog/references/staged.tsv",
        decoys="data/results/phase4_catalog/references/decoy_candidates.tsv",
        code=["scripts/prepare_expanded_references.py", "scripts/assembly_completion.py"],
    output:
        references=f"{NEW}/nominated.tsv",
        decoys=f"{NEW}/decoy_candidates.tsv",
        fixed_decoys=f"{REF}/decoy_candidates.tsv",
    resources: mem_mb=2000, runtime=60, partition="ycga"
    shell: "{PY:q} scripts/prepare_expanded_references.py select"

rule collect_new_references:
    input:
        references=rules.select_new_references.output.references,
        decoys=rules.select_new_references.output.decoys,
        archive="data/results/phase4_catalog/oceandna/fasta_species-representatives.tar",
        code="workflow/scripts/collect_reference_genomes.py",
    output:
        genomes=directory(f"{NEW}/genomes"),
        table=f"{NEW}/staged.tsv",
        batch=f"{NEW}/oceandna_batch.tsv",
    params: genomes=P4["references"]["gtdb_genomes"]
    resources: mem_mb=4000, runtime=120, partition="ycga"
    shell:
        "{PY:q} {input.code:q} --references {input.references:q} --decoys {input.decoys:q} --gtdb-genomes {params.genomes:q} --oceandna-archive {input.archive:q} --genomes {output.genomes:q} --table {output.table:q} --oceandna-batch {output.batch:q}"

rule check_new_references_ncbi:
    input:
        staged=rules.collect_new_references.output.table,
        code="workflow/scripts/check_ncbi_status.py",
    output: f"{NEW}/ncbi_status.tsv"
    resources: mem_mb=2000, runtime=120, partition="ycga"
    shell: "{PY:q} {input.code:q} --genomes {input.staged:q} --email casey.dunn@yale.edu --output {output:q}"

rule fetch_new_substitutes:
    input:
        status=rules.check_new_references_ncbi.output,
        code="workflow/scripts/fetch_ncbi_substitutes.py",
    output:
        genomes=directory(f"{NEW}/substitutes"),
        table=f"{NEW}/substitutes.tsv",
    resources: mem_mb=2000, runtime=120, partition="ycga"
    shell: "{PY:q} {input.code:q} --status {input.status:q} --genomes {output.genomes:q} --table {output.table:q}"

rule assess_new_references_checkm2:
    input:
        genomes=rules.collect_new_references.output.genomes,
        substitutes=rules.fetch_new_substitutes.output.genomes,
        database=P4["databases"]["checkm2"],
    output: f"{NEW}/checkm2_quality.tsv"
    params: tools=P4["tool_prefixes"]["checkm2"] + "/bin"
    threads: 16
    resources: mem_mb=64000, runtime=480, partition="ycga"
    log: f"{WORK}/logs/expanded_references/checkm2.log"
    shell:
        r"""
        export PATH={params.tools:q}:$PATH
        tmp=$(mktemp -d "${{TMPDIR:-/tmp}}/expanded-ref-qc.XXXXXX")
        trap 'rm -rf "$tmp"' EXIT
        mkdir "$tmp/genomes"
        for f in {input.genomes:q}/*.fa {input.substitutes:q}/*.fa; do
          if [ -e "$f" ]; then ln -s "$(realpath "$f")" "$tmp/genomes/"; fi
        done
        checkm2 predict --threads {threads} --input "$tmp/genomes" --extension fa \
          --database_path {input.database:q} --output-directory "$tmp/out" --force > {log:q} 2>&1
        cp "$tmp/out/quality_report.tsv" {output:q}
        """

rule classify_new_oceandna_gtdbtk:
    input:
        batch=rules.collect_new_references.output.batch,
        genomes=rules.collect_new_references.output.genomes,
        reference=ancient("/vast/palmer/scratch/dunn/cwd7/siph_associates/phase3_catalog/gtdbtk_r220"),
        mash=ancient("/vast/palmer/scratch/dunn/cwd7/siph_associates/phase3_catalog/gtdbtk_mash/gtdb_r220.msh"),
    output:
        directory(f"{NEW}/gtdb"),
    params: tools=P4["tool_prefixes"]["gtdbtk"] + "/bin"
    threads: 16
    resources: mem_mb=150000, runtime=720, partition="ycga"
    log: f"{WORK}/logs/expanded_references/gtdb.log"
    shell:
        "export PATH={params.tools:q}:$PATH GTDBTK_DATA_PATH={input.reference:q}; gtdbtk classify_wf --batchfile {input.batch:q} --out_dir {output:q} --cpus {threads} --pplacer_cpus 8 --mash_db {input.mash:q} > {log:q} 2>&1"

rule new_reference_evidence:
    input:
        rules.assess_new_references_checkm2.output,
        rules.classify_new_oceandna_gtdbtk.output,
    default_target: True
