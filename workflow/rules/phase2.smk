rule validate_screen_sample:
    input:
        trim_r1=f"{SCRATCH}/trimmed/{{sample}}_R1.fastq.gz",
        trim_r2=f"{SCRATCH}/trimmed/{{sample}}_R2.fastq.gz",
        fastp=f"{WORK}/trim/{{sample}}.fastp.json",
        trim_provenance=f"{WORK}/provenance/trim/{{sample}}.json",
        kraken=f"{WORK}/screens/kraken/{{sample}}.report.tsv",
        genus=f"{WORK}/screens/bracken/{{sample}}.G.tsv",
        species=f"{WORK}/screens/bracken/{{sample}}.S.tsv",
        sylph=f"{WORK}/screens/sylph/{{sample}}.profile.tsv",
        phyloflash=f"{WORK}/screens/phyloflash/{{sample}}.tar.gz",
        kraken_provenance=f"{WORK}/provenance/kraken_bracken/{{sample}}.json",
        sylph_provenance=f"{WORK}/provenance/sylph/{{sample}}.json",
        phyloflash_provenance=f"{WORK}/provenance/phyloflash/{{sample}}.json",
    output:
        f"{WORK}/validation/screens/{{sample}}.json"
    params:
        library_id=lambda wildcards: SAMPLES[wildcards.sample]["library_id"],
        specimen_id=lambda wildcards: SAMPLES[wildcards.sample]["specimen_id"],
        study=lambda wildcards: SAMPLES[wildcards.sample]["study"],
        route=host_route,
        raw_pairs=lambda wildcards: int(SAMPLES[wildcards.sample]["read_pairs"]),
        cap_pairs=int(config["parameters"]["cap_pairs"]),
    log:
        f"{WORK}/logs/validation/screens/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/validation/screens/{{sample}}.tsv",
    threads: 1
    priority: 200
    resources:
        mem_mb=config["resources"]["validation"]["mem_mb"],
        runtime=config["resources"]["validation"]["runtime"],
        partition=config["resources"]["validation"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/validate_screen_sample.py"


rule aggregate_phase2:
    input:
        validations=lambda wildcards: phase2_validation_paths(phase2_scope_samples(wildcards)),
        genus=lambda wildcards: expand(f"{WORK}/screens/bracken/{{sample}}.G.tsv", sample=phase2_scope_samples(wildcards)),
        species=lambda wildcards: expand(f"{WORK}/screens/bracken/{{sample}}.S.tsv", sample=phase2_scope_samples(wildcards)),
        sylph=lambda wildcards: expand(f"{WORK}/screens/sylph/{{sample}}.profile.tsv", sample=phase2_scope_samples(wildcards)),
        phyloflash=lambda wildcards: expand(f"{WORK}/screens/phyloflash/{{sample}}.tar.gz", sample=phase2_scope_samples(wildcards)),
    output:
        qc=f"{WORK}/aggregation/{{scope}}/library_qc.tsv",
        nominations=f"{WORK}/aggregation/{{scope}}/candidate_nominations.tsv",
        provenance=f"{WORK}/provenance/aggregation/{{scope}}.json",
    params:
        sample_ids=phase2_scope_samples,
    wildcard_constraints:
        scope="pilot|cohort"
    log:
        f"{WORK}/logs/aggregation/{{scope}}.log",
    benchmark:
        f"{WORK}/benchmarks/aggregation/{{scope}}.tsv",
    threads: 1
    priority: 300
    resources:
        mem_mb=config["resources"]["aggregation"]["mem_mb"],
        runtime=config["resources"]["aggregation"]["runtime"],
        partition=config["resources"]["aggregation"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/aggregate_phase2.py"


rule screen_pilot:
    input:
        validations=phase2_validation_paths(PILOT_IDS),
        qc=f"{WORK}/aggregation/pilot/library_qc.tsv",
        nominations=f"{WORK}/aggregation/pilot/candidate_nominations.tsv",
        provenance=f"{WORK}/provenance/aggregation/pilot.json",
    output:
        f"{WORK}/stages/screen_pilot.done"
    params:
        expected_samples=len(PILOT_IDS),
        expected_studies=len({SAMPLES[s]["study"] for s in PILOT_IDS}),
        expected_routes=len({SAMPLES[s]["host_route"] for s in PILOT_IDS}),
    log:
        f"{WORK}/logs/validation/screen_pilot.log",
    benchmark:
        f"{WORK}/benchmarks/validation/screen_pilot.tsv",
    threads: 1
    resources:
        mem_mb=config["resources"]["validation"]["mem_mb"],
        runtime=config["resources"]["validation"]["runtime"],
        partition=config["resources"]["validation"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/validate_phase2_stage.py"


rule screen_cohort:
    input:
        validations=phase2_validation_paths(SAMPLE_IDS),
        qc=f"{WORK}/aggregation/cohort/library_qc.tsv",
        nominations=f"{WORK}/aggregation/cohort/candidate_nominations.tsv",
        provenance=f"{WORK}/provenance/aggregation/cohort.json",
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
        mem_mb=config["resources"]["validation"]["mem_mb"],
        runtime=config["resources"]["validation"]["runtime"],
        partition=config["resources"]["validation"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/validate_phase2_stage.py"
