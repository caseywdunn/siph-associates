rule snapshot_run_configuration:
    input:
        config=lambda wildcards: config["config_path"],
        samples=config["samples"],
        manifest=config["manifest"],
        freeze=config["manifest_freeze"],
    output:
        config=f"{WORK}/provenance/run/config.yaml",
        samples=f"{WORK}/provenance/run/samples.tsv",
        manifest_freeze=f"{WORK}/provenance/run/manifest.freeze.sha256",
        checksums=f"{WORK}/provenance/run/inputs.sha256",
    log:
        f"{WORK}/logs/validation/snapshot_run.log",
    benchmark:
        f"{WORK}/benchmarks/validation/snapshot_run.tsv",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30,
        partition="day",
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/snapshot_run.py"


rule phase1_ready:
    input:
        samples=config["samples"],
        manifest=config["manifest"],
        freeze=config["manifest_freeze"],
        workflow="Snakefile",
        config="config/config.yaml",
        snapshot=f"{WORK}/provenance/run/inputs.sha256",
    output:
        f"{WORK}/stages/phase1.ready"
    params:
        expected_samples=len(SAMPLE_IDS),
    log:
        f"{WORK}/logs/validation/phase1_ready.log",
    benchmark:
        f"{WORK}/benchmarks/validation/phase1_ready.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
        partition="day",
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/validate_workflow_inputs.py"


rule validate_sample:
    input:
        trim_r1=f"{SCRATCH}/trimmed/{{sample}}_R1.fastq.gz",
        trim_r2=f"{SCRATCH}/trimmed/{{sample}}_R2.fastq.gz",
        fastp=f"{WORK}/trim/{{sample}}.fastp.json",
        kraken=f"{WORK}/screens/kraken/{{sample}}.report.tsv",
        genus=f"{WORK}/screens/bracken/{{sample}}.G.tsv",
        species=f"{WORK}/screens/bracken/{{sample}}.S.tsv",
        sylph=f"{WORK}/screens/sylph/{{sample}}.profile.tsv",
        phyloflash=f"{WORK}/screens/phyloflash/{{sample}}.tar.gz",
        host_r1=f"{WORK}/host_handling/{{sample}}_R1.fastq.gz",
        host_r2=f"{WORK}/host_handling/{{sample}}_R2.fastq.gz",
        host_metrics=f"{WORK}/host_handling/{{sample}}.metrics.json",
    output:
        f"{WORK}/validation/samples/{{sample}}.json"
    params:
        study=lambda wildcards: SAMPLES[wildcards.sample]["study"],
        route=host_route,
    log:
        f"{WORK}/logs/validation/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/validation/{{sample}}.tsv",
    threads: 1
    resources:
        mem_mb=config["resources"]["validation"]["mem_mb"],
        runtime=config["resources"]["validation"]["runtime"],
        partition=config["resources"]["validation"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/validate_sample.py"


rule phase1_smoke:
    input:
        expand(f"{WORK}/validation/samples/{{sample}}.json", sample=SAMPLE_IDS)
    output:
        f"{WORK}/stages/phase1_smoke.done"
    params:
        expected_samples=len(SAMPLE_IDS),
        expected_studies=len({SAMPLES[s]["study"] for s in SAMPLE_IDS}),
        expected_routes=len({SAMPLES[s]["host_route"] for s in SAMPLE_IDS}),
    log:
        f"{WORK}/logs/validation/phase1_smoke.log",
    benchmark:
        f"{WORK}/benchmarks/validation/phase1_smoke.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
        partition="day",
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/validate_stage.py"


rule screen_cohort:
    input:
        expand(f"{WORK}/validation/samples/{{sample}}.json", sample=SAMPLE_IDS)
    output:
        f"{WORK}/stages/screen_cohort.done"
    params:
        expected_samples=len(SAMPLE_IDS),
        expected_studies=len({SAMPLES[s]["study"] for s in SAMPLE_IDS}),
        expected_routes=len({SAMPLES[s]["host_route"] for s in SAMPLE_IDS}),
    log:
        f"{WORK}/logs/validation/screen_cohort.log",
    benchmark:
        f"{WORK}/benchmarks/validation/screen_cohort.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
        partition="day",
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/validate_stage.py"
