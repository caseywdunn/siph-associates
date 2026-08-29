rule make_fixture_fastq:
    input:
        r1=source_r1,
        r2=source_r2,
    output:
        r1=f"{FIXTURE_ROOT}/raw/{{sample}}_R1.fastq.gz",
        r2=f"{FIXTURE_ROOT}/raw/{{sample}}_R2.fastq.gz",
        metadata=f"{FIXTURE_ROOT}/raw/{{sample}}.fixture.json",
    params:
        pairs=lambda wildcards: int(config["fixture_pairs"]),
        enabled=lambda wildcards: bool(config["fixture_enabled"]),
    log:
        f"{WORK}/logs/make_fixture/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/make_fixture/{{sample}}.tsv",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30,
        partition="day",
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/make_fixture.py"


rule trim_reads:
    input:
        r1=raw_r1,
        r2=raw_r2,
        run_snapshot=f"{WORK}/provenance/run/inputs.sha256",
    output:
        r1=temp(f"{SCRATCH}/trimmed/{{sample}}_R1.fastq.gz") if config.get("cleanup_trimmed_after_validation") else f"{SCRATCH}/trimmed/{{sample}}_R1.fastq.gz",
        r2=temp(f"{SCRATCH}/trimmed/{{sample}}_R2.fastq.gz") if config.get("cleanup_trimmed_after_validation") else f"{SCRATCH}/trimmed/{{sample}}_R2.fastq.gz",
        fastp_json=f"{WORK}/trim/{{sample}}.fastp.json",
        fastp_html=f"{WORK}/trim/{{sample}}.fastp.html",
        provenance=f"{WORK}/provenance/trim/{{sample}}.json",
    params:
        expected_pairs=lambda wildcards: (
            int(config["fixture_pairs"]) if config["fixture_enabled"]
            else int(SAMPLES[wildcards.sample]["read_pairs"])
        ),
        cap_pairs=lambda wildcards: int(config["parameters"]["cap_pairs"]),
        cap_seed=lambda wildcards: int(config["parameters"]["cap_seed"]),
    log:
        f"{WORK}/logs/trim/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/trim/{{sample}}.tsv",
    threads: config["resources"]["trim"]["threads"]
    priority: 0
    resources:
        mem_mb=config["resources"]["trim"]["mem_mb"],
        runtime=config["resources"]["trim"]["runtime"],
        partition=config["resources"]["trim"]["partition"],
        cap_slots=lambda wildcards: (
            1 if int(SAMPLES[wildcards.sample]["read_pairs"]) > int(config["parameters"]["cap_pairs"])
            else 0
        ),
    envmodules:
        "fastp/0.23.2-GCCcore-10.2.0",
        "BBMap/38.90-GCCcore-10.2.0"
    conda:
        "../../envs/fastp.yaml"
    script:
        "../scripts/trim_reads.py"
