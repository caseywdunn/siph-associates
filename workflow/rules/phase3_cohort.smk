rule phase3_prepare_cohort_input:
    input:
        r1=ancient(f"{SCRATCH}/trimmed/{{sample}}_R1.fastq.gz"),
        r2=ancient(f"{SCRATCH}/trimmed/{{sample}}_R2.fastq.gz"),
        upstream=f"{WORK}/provenance/trim/{{sample}}.json",
        index=host_index_inputs,
    output:
        r1=f"{WORK}/phase3_cohort/inputs/{{sample}}_R1.fastq.gz",
        r2=f"{WORK}/phase3_cohort/inputs/{{sample}}_R2.fastq.gz",
        metrics=f"{WORK}/phase3_cohort/inputs/{{sample}}.json",
        provenance=f"{WORK}/provenance/phase3_cohort/inputs/{{sample}}.json",
    params:
        sample=lambda wildcards: wildcards.sample,
        route=host_route,
        strategy=phase3_cohort_strategy,
        prefix=config["tool_prefixes"]["assembly"],
        index_prefix=lambda wildcards: f"{WORK}/reference_indices/{host_route(wildcards)}/host",
        fixed_pairs=int(PHASE3_COHORT["fixed_effort_pairs"]),
        fixed_seed=int(PHASE3_COHORT["fixed_effort_seed"]),
    log:
        f"{WORK}/logs/phase3_cohort/inputs/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_cohort/inputs/{{sample}}.tsv",
    threads: PHASE3_COHORT["resources"]["input_preparation"]["threads"]
    resources:
        mem_mb=PHASE3_COHORT["resources"]["input_preparation"]["mem_mb"],
        runtime=PHASE3_COHORT["resources"]["input_preparation"]["runtime"],
        partition=PHASE3_COHORT["resources"]["input_preparation"]["partition"],
    conda:
        "../../envs/fastp.yaml"
    script:
        "../scripts/phase3_prepare_cohort_input.py"


rule aggregate_phase3_inputs:
    input:
        config=str(PHASE3_COHORT_CONFIG_PATH),
        metrics=expand(f"{WORK}/phase3_cohort/inputs/{{sample}}.json", sample=SAMPLE_IDS),
    output:
        table=f"{WORK}/phase3_cohort/input_eligibility.tsv",
        provenance=f"{WORK}/provenance/phase3_cohort/input_eligibility.json",
    params:
        sample_ids=SAMPLE_IDS,
        minimum_pairs=int(PHASE3_COHORT["minimum_assembly_input_pairs"]),
    log:
        f"{WORK}/logs/phase3_cohort/input_eligibility.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_cohort/input_eligibility.tsv",
    threads: 1
    resources:
        mem_mb=PHASE3_COHORT["resources"]["validation"]["mem_mb"],
        runtime=PHASE3_COHORT["resources"]["validation"]["runtime"],
        partition=PHASE3_COHORT["resources"]["validation"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/aggregate_phase3_inputs.py"


rule phase3_inputs:
    input:
        table=f"{WORK}/phase3_cohort/input_eligibility.tsv",
        provenance=f"{WORK}/provenance/phase3_cohort/input_eligibility.json",
    output:
        f"{WORK}/stages/phase3_inputs.done"
    params:
        expected_samples=len(SAMPLE_IDS),
        expected_routes=len({SAMPLES[s]["host_route"] for s in SAMPLE_IDS}),
    log:
        f"{WORK}/logs/phase3_cohort/validate_inputs.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_cohort/validate_inputs.tsv",
    threads: 1
    resources:
        mem_mb=PHASE3_COHORT["resources"]["validation"]["mem_mb"],
        runtime=PHASE3_COHORT["resources"]["validation"]["runtime"],
        partition=PHASE3_COHORT["resources"]["validation"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/validate_phase3_inputs.py"
