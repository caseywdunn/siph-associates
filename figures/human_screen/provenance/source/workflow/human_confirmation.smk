"""Search representative stringent human-compatible candidates against broad nt.

Only completed human-screen outputs are consumed; this entry point cannot rerun
raw-read screening or any existing microbial/viral mapping analysis.
"""
import csv
import json
from pathlib import Path

with open("config/human_screen.json") as handle:
    HUMAN = json.load(handle)
SCOPE = config.get("scope", "pilot")
if SCOPE not in ("pilot", "cohort"):
    raise ValueError("scope must be pilot or cohort")
ROOT = HUMAN["run"]
OUT = ROOT + "/confirmation/" + SCOPE
PYTHON = HUMAN["python"]
if SCOPE == "pilot":
    SAMPLES = HUMAN["pilot_samples"]
else:
    with open(HUMAN["manifest"]) as handle:
        SAMPLES = [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")]

rule prepare_human_candidate_queries:
    input:
        fasta=expand(ROOT + "/samples/{sample}/candidates.fasta", sample=SAMPLES),
        pairs=expand(ROOT + "/samples/{sample}/pairs.tsv", sample=SAMPLES),
        config="config/human_screen.json",
        script="scripts/human_screen/prepare_confirmation.py",
    output:
        fasta=OUT + "/queries.fasta",
        sources=OUT + "/query_sources.tsv",
        selection=OUT + "/selection.json",
    shell:
        "{PYTHON} {input.script:q} --config {input.config:q} --scope {SCOPE} --max-pairs 500 --output {OUT:q}"

rule search_human_candidates_blastn:
    input:
        fasta=rules.prepare_human_candidate_queries.output.fasta,
        database="/gpfs/gibbs/data/db/blast/nt.nal",
    output: OUT + "/nt_hits.tsv"
    threads: 16
    resources: mem_mb=60000
    log: OUT + "/blastn.log"
    shell:
        """
        if [ -s {input.fasta:q} ]; then
          blastn -task megablast -query {input.fasta:q} -db /gpfs/gibbs/data/db/blast/nt \\
            -perc_identity 98 -qcov_hsp_perc 95 -max_target_seqs 100 -max_hsps 1 \\
            -evalue 1e-15 -num_threads {threads} \\
            -outfmt '6 qseqid saccver pident length qlen qstart qend evalue bitscore staxids sscinames' \\
            -out {output:q} > {log:q} 2>&1
        else
          touch {output:q} {log:q}
        fi
        """

rule review_human_candidate_taxonomy:
    input:
        hits=rules.search_human_candidates_blastn.output,
        sources=rules.prepare_human_candidate_queries.output.sources,
        script="scripts/human_screen/review_nt.py",
        names="/gpfs/gibbs/data/db/blast/names.dmp",
    output:
        table=OUT + "/taxonomic_review.tsv",
        summary=OUT + "/taxonomic_review.json",
        pairs=OUT + "/pair_review.tsv",
    shell: "{PYTHON} {input.script:q} --hits {input.hits:q} --sources {input.sources:q} --output {output.table:q} --names {input.names:q}"

rule stage_human_screen_report:
    input:
        samples=ROOT + "/cohort_summary.tsv",
        controls=ROOT + "/controls.tsv",
        control_check=ROOT + "/controls.passed.json",
        nt=rules.review_human_candidate_taxonomy.output.summary,
        selection=rules.prepare_human_candidate_queries.output.selection,
        script="scripts/human_screen/stage_report.py",
        config="config/human_screen.json",
    output:
        samples="figures/human_screen/sample_summary.tsv",
        controls="figures/human_screen/control_summary.tsv",
        nt="figures/human_screen/nt_summary.json",
        summary="figures/human_screen/summary.json",
        values="figures/human_screen/values.tex",
        provenance="figures/human_screen/provenance.json",
    params:
        executions=lambda w: config.get("execution_ids", "").split(","),
    shell:
        "{PYTHON} {input.script:q} --config {input.config:q} --confirmation-scope {SCOPE} --execution-ids {params.executions:q} --output figures/human_screen"

rule report:
    input: rules.stage_human_screen_report.output

rule all:
    input: rules.review_human_candidate_taxonomy.output
    default_target: True
