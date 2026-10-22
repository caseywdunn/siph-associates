"""Assess human-compatible read pairs without altering associate catalogs.

Per-source-file prefixes and known-origin controls -> QC -> stringent mapping
-> explicit postfilters -> pilot/cohort summaries. Cohort execution requires
successful controls. Candidate sequence taxonomy remains a separate review.
"""
import csv
from pathlib import Path

configfile: "config/human_screen.json"
RUN = config["run"]
CONFIG = "config/human_screen.json"
PYTHON = config["python"]
with open(config["manifest"]) as handle:
    ROWS = {r["sample_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
SAMPLES = list(ROWS)
CONTROLS = ["control_human"] + ["control_" + key for key in config["negative_references"]]

rule prepare_human_screen_reads:
    input:
        script="scripts/human_screen/prepare_reads.py",
        config=CONFIG,
        manifest=config["manifest"],
        controls=lambda w: [] if w.sample.startswith("control_") else [RUN + "/controls.passed.json"],
        source=lambda w: ([config["human_fasta"] if w.sample == "control_human" else config["negative_references"][w.sample.removeprefix("control_")]] if w.sample.startswith("control_") else ROWS[w.sample]["r1_paths"].split(";") + ROWS[w.sample]["r2_paths"].split(";")),
    output:
        r1=temp(RUN + "/reads/{sample}_R1.fastq.gz"),
        r2=temp(RUN + "/reads/{sample}_R2.fastq.gz"),
        metrics=RUN + "/reads/{sample}.json",
    threads: 1
    resources: mem_mb=2000
    shell:
        "{PYTHON} {input.script:q} --config {input.config:q} --sample {wildcards.sample:q} --r1 {output.r1:q} --r2 {output.r2:q} --metrics {output.metrics:q}"

rule trim_human_screen_fastp:
    input:
        r1=rules.prepare_human_screen_reads.output.r1,
        r2=rules.prepare_human_screen_reads.output.r2,
    output:
        r1=temp(RUN + "/trim/{sample}_R1.fastq.gz"),
        r2=temp(RUN + "/trim/{sample}_R2.fastq.gz"),
        json=RUN + "/trim/{sample}.json",
        html=RUN + "/trim/{sample}.html",
    threads: 2
    resources: mem_mb=2000
    log: RUN + "/logs/{sample}.fastp.log"
    shell:
        "fastp -i {input.r1:q} -I {input.r2:q} -o {output.r1:q} -O {output.r2:q} --detect_adapter_for_pe --dont_eval_duplication --length_required 100 --thread {threads} --json {output.json:q} --html {output.html:q} > {log:q} 2>&1"

rule map_human_screen_bowtie2:
    input:
        r1=rules.trim_human_screen_fastp.output.r1,
        r2=rules.trim_human_screen_fastp.output.r2,
        index=[config["human_index"] + suffix for suffix in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")],
    output: bam=RUN + "/mapping/{sample}.bam"
    params: index=config["human_index"]
    threads: 4
    resources: mem_mb=7000
    benchmark: RUN + "/benchmarks/{sample}.mapping.tsv"
    log: RUN + "/logs/{sample}.bowtie2.log"
    shell:
        """
        bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --no-unal --reorder \\
          --score-min L,0,-0.6 --seed 20261022 -I 0 -X 1000 -p {threads} \\
          -x {params.index:q} -1 {input.r1:q} -2 {input.r2:q} 2> {log:q} | samtools view -b -o {output.bam:q} -
        """

rule filter_human_screen_alignments:
    input:
        bam=rules.map_human_screen_bowtie2.output.bam,
        selected=rules.prepare_human_screen_reads.output.metrics,
        fastp=rules.trim_human_screen_fastp.output.json,
        script="scripts/human_screen/filter_alignments.py",
        config=CONFIG,
    output:
        summary=RUN + "/samples/{sample}/summary.json",
        pairs=RUN + "/samples/{sample}/pairs.tsv",
        fasta=RUN + "/samples/{sample}/candidates.fasta",
        r1=RUN + "/samples/{sample}/candidates_R1.fastq.gz",
        r2=RUN + "/samples/{sample}/candidates_R2.fastq.gz",
    params: output=lambda w: RUN + "/samples/" + w.sample
    threads: 1
    resources: mem_mb=2000
    shell:
        "{PYTHON} {input.script:q} --config {input.config:q} --sample {wildcards.sample:q} --bam {input.bam:q} --input-metrics {input.selected:q} --fastp-metrics {input.fastp:q} --output {params.output:q}"

rule validate_human_screen_controls:
    input:
        summaries=expand(RUN + "/samples/{sample}/summary.json", sample=CONTROLS),
        script="scripts/human_screen/report.py",
    output:
        table=RUN + "/controls.tsv",
        passed=RUN + "/controls.passed.json",
    shell: "{PYTHON} {input.script:q} --controls --output {output.table:q} {input.summaries:q}"

rule report_human_screen_pilot:
    input:
        summaries=expand(RUN + "/samples/{sample}/summary.json", sample=config["pilot_samples"]),
        controls=rules.validate_human_screen_controls.output.passed,
        script="scripts/human_screen/report.py",
    output: RUN + "/pilot_summary.tsv"
    shell: "{PYTHON} {input.script:q} --output {output:q} {input.summaries:q}"

rule report_human_screen_cohort:
    input:
        summaries=expand(RUN + "/samples/{sample}/summary.json", sample=SAMPLES),
        controls=rules.validate_human_screen_controls.output.passed,
        script="scripts/human_screen/report.py",
    output: RUN + "/cohort_summary.tsv"
    shell: "{PYTHON} {input.script:q} --output {output:q} {input.summaries:q}"

rule pilot:
    input: rules.report_human_screen_pilot.output

rule controls:
    input: rules.validate_human_screen_controls.output

rule all:
    input: rules.report_human_screen_cohort.output
    default_target: True
