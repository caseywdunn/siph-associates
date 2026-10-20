"""Screen all remaining libraries using capped quality-filtered whole reads."""

import csv

configfile: "config/fish_12s_full_reads.json"

ROOT, PY = config["run"], config["python"]
READ_LABEL = config["read_pair_id"]
with open(config["libraries"]) as handle:
    LIBRARIES = list(csv.DictReader(handle, delimiter="\t"))
SAMPLES = {r["analysis_sample_id"]: r for r in LIBRARIES}
if len(SAMPLES) != len(LIBRARIES):
    raise ValueError("Duplicate fish-screen library identifiers")

wildcard_constraints:
    sample="[^/]+",


def fish_read(w, mate):
    return f"{ROOT}/trimmed/{w.sample}_R{mate}.fastq.gz"


def fish_manifest(w):
    return f"{ROOT}/inputs/{w.sample}.tsv"


def fish_index(w):
    return 0


rule prepare_fish_full_reads:
    input:
        r1=lambda w: SAMPLES[w.sample]["r1_paths"].split(";"),
        r2=lambda w: SAMPLES[w.sample]["r2_paths"].split(";"),
        run_snapshot=ancient("data/results/provenance/run/inputs.sha256"),
    output:
        r1=temp(f"{ROOT}/trimmed/{{sample}}_R1.fastq.gz"),
        r2=temp(f"{ROOT}/trimmed/{{sample}}_R2.fastq.gz"),
        fastp_json=f"{ROOT}/trim/{{sample}}.fastp.json",
        fastp_html=f"{ROOT}/trim/{{sample}}.fastp.html",
        provenance=f"{ROOT}/trim/{{sample}}.json",
    params:
        expected_pairs=lambda w: int(SAMPLES[w.sample]["raw_pairs"]),
        cap_pairs=int(config["cap_pairs"]),
        cap_seed=int(config["cap_seed"]),
    threads: 8
    resources:
        mem_mb=32000,
    log:
        f"{ROOT}/logs/trim/{{sample}}.log",
    benchmark:
        f"{ROOT}/benchmarks/trim/{{sample}}.tsv",
    conda:
        "../envs/fastp.yaml"
    script:
        "scripts/trim_reads.py"


rule inventory_fish_full_reads:
    input:
        r1=rules.prepare_fish_full_reads.output.r1,
        r2=rules.prepare_fish_full_reads.output.r2,
        provenance=rules.prepare_fish_full_reads.output.provenance,
        libraries=config["libraries"],
        code="scripts/fish_12s/full_read_inputs.py",
    output:
        f"{ROOT}/inputs/{{sample}}.tsv",
    shell:
        "{PY:q} {input.code:q} --libraries {input.libraries:q} --sample {wildcards.sample:q} "
        "--provenance {input.provenance:q} --output {output:q}"


include: "rules/fish_12s_search.smk"


rule review_fish_full_reads:
    input:
        hits=f"{ROOT}/union/{{sample}}/{READ_LABEL}/alignment/hits.tsv",
        fasta=f"{ROOT}/union/{{sample}}/{READ_LABEL}/candidates.fasta",
        pairs=f"{ROOT}/union/{{sample}}/{READ_LABEL}/candidate_pairs.tsv",
        run=f"{ROOT}/union/{{sample}}/{READ_LABEL}/run.json",
        manifest=fish_manifest,
        metadata=config["reference_metadata"],
        code=["scripts/fish_12s/summarize_cohort.py", "scripts/fish_12s/classify_candidates.py"],
    output:
        summary=f"{ROOT}/review/{{sample}}/sample_summary.tsv",
        pairs=f"{ROOT}/review/{{sample}}/pair_review.tsv",
        run=f"{ROOT}/review/{{sample}}/run.json",
    params:
        config=lambda w: config.get("config_path", "config/fish_12s_full_reads.json"),
        directory=lambda w: f"{ROOT}/review/{w.sample}",
    shell:
        "{PY:q} scripts/fish_12s/summarize_cohort.py --config {params.config:q} "
        "--manifest {input.manifest:q} --reference-metadata {input.metadata:q} --output {params.directory:q}"


rule full_reads:
    input:
        expand(f"{ROOT}/review/{{sample}}/sample_summary.tsv", sample=sorted(SAMPLES)),
    default_target: True
