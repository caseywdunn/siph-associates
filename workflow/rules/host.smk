rule make_fixture_reference:
    input:
        lambda wildcards: config["fixture_reference_sources"][wildcards.route]
    output:
        f"{FIXTURE_ROOT}/references/{{route}}.fasta"
    params:
        enabled=lambda wildcards: bool(config["fixture_enabled"]),
        bases=lambda wildcards: int(config.get("fixture_reference_bases", 1000000)),
    log:
        f"{WORK}/logs/fixture_reference/{{route}}.log",
    benchmark:
        f"{WORK}/benchmarks/fixture_reference/{{route}}.tsv",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30,
        partition="day",
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/make_fixture_reference.py"


rule index_host_reference:
    input:
        reference=reference_source,
        run_snapshot=f"{WORK}/provenance/run/inputs.sha256",
    output:
        amb=f"{WORK}/reference_indices/{{route}}/host.amb",
        ann=f"{WORK}/reference_indices/{{route}}/host.ann",
        bwt=f"{WORK}/reference_indices/{{route}}/host.bwt",
        pac=f"{WORK}/reference_indices/{{route}}/host.pac",
        sa=f"{WORK}/reference_indices/{{route}}/host.sa",
        provenance=f"{WORK}/provenance/reference_index/{{route}}.json",
    params:
        prefix=config["tool_prefixes"]["assembly"],
        index_prefix=lambda wildcards: f"{WORK}/reference_indices/{wildcards.route}/host",
    log:
        f"{WORK}/logs/reference_index/{{route}}.log",
    benchmark:
        f"{WORK}/benchmarks/reference_index/{{route}}.tsv",
    threads: 1
    resources:
        mem_mb=64000,
        runtime=1440,
        partition="day",
    conda:
        "../../envs/host-mapping.yaml"
    script:
        "../scripts/index_reference.py"


rule host_handling:
    input:
        r1=ancient(f"{SCRATCH}/trimmed/{{sample}}_R1.fastq.gz"),
        r2=ancient(f"{SCRATCH}/trimmed/{{sample}}_R2.fastq.gz"),
        classified_r1=f"{WORK}/screens/kraken/{{sample}}.classified_1.fastq.gz",
        classified_r2=f"{WORK}/screens/kraken/{{sample}}.classified_2.fastq.gz",
        trim_provenance=f"{WORK}/provenance/trim/{{sample}}.json",
        kraken_provenance=f"{WORK}/provenance/kraken_bracken/{{sample}}.json",
        index=host_index_inputs,
    output:
        r1=f"{WORK}/host_handling/{{sample}}_R1.fastq.gz",
        r2=f"{WORK}/host_handling/{{sample}}_R2.fastq.gz",
        metrics=f"{WORK}/host_handling/{{sample}}.metrics.json",
        provenance=f"{WORK}/provenance/host_handling/{{sample}}.json",
    params:
        route=host_route,
        prefix=config["tool_prefixes"]["assembly"],
        index_prefix=lambda wildcards: f"{WORK}/reference_indices/{host_route(wildcards)}/host",
    log:
        f"{WORK}/logs/host_handling/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/host_handling/{{sample}}.tsv",
    threads: config["resources"]["host_handling"]["threads"]
    resources:
        mem_mb=config["resources"]["host_handling"]["mem_mb"],
        runtime=config["resources"]["host_handling"]["runtime"],
        partition=config["resources"]["host_handling"]["partition"],
    conda:
        "../../envs/host-mapping.yaml"
    script:
        "../scripts/host_handling.py"
