rule benchmark_trim_reads:
    input:
        r1=raw_r1,
        r2=raw_r2,
        run_snapshot=f"{WORK}/stages/screen_cohort.done",
    output:
        r1=temp(f"{SCRATCH}/phase3_benchmark/trimmed/{{sample}}_R1.fastq.gz"),
        r2=temp(f"{SCRATCH}/phase3_benchmark/trimmed/{{sample}}_R2.fastq.gz"),
        fastp_json=f"{WORK}/phase3_benchmark/trim/{{sample}}.fastp.json",
        fastp_html=f"{WORK}/phase3_benchmark/trim/{{sample}}.fastp.html",
        provenance=f"{WORK}/provenance/phase3_benchmark/trim/{{sample}}.json",
    params:
        expected_pairs=lambda wildcards: int(SAMPLES[wildcards.sample]["read_pairs"]),
        cap_pairs=int(config["parameters"]["cap_pairs"]),
        cap_seed=int(config["parameters"]["cap_seed"]),
    log:
        f"{WORK}/logs/phase3_benchmark/trim/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_benchmark/trim/{{sample}}.tsv",
    threads: config["resources"]["trim"]["threads"]
    resources:
        mem_mb=48000,
        runtime=360,
        partition="day",
        cap_slots=lambda wildcards: 1 if int(SAMPLES[wildcards.sample]["read_pairs"]) > int(config["parameters"]["cap_pairs"]) else 0,
    envmodules:
        "fastp/0.23.2-GCCcore-10.2.0",
        "BBMap/38.90-GCCcore-10.2.0"
    conda:
        "../../envs/fastp.yaml"
    script:
        "../scripts/trim_reads.py"


rule benchmark_fixed_effort_reads:
    input:
        r1=lambda wildcards: f"{SCRATCH}/phase3_benchmark/trimmed/{phase3_sample(wildcards)}_R1.fastq.gz",
        r2=lambda wildcards: f"{SCRATCH}/phase3_benchmark/trimmed/{phase3_sample(wildcards)}_R2.fastq.gz",
        upstream=lambda wildcards: f"{WORK}/provenance/phase3_benchmark/trim/{phase3_sample(wildcards)}.json",
    output:
        r1=temp(f"{SCRATCH}/phase3_benchmark/inputs/fixed/{{benchmark_id}}_R1.fastq.gz"),
        r2=temp(f"{SCRATCH}/phase3_benchmark/inputs/fixed/{{benchmark_id}}_R2.fastq.gz"),
        provenance=f"{WORK}/provenance/phase3_benchmark/inputs/fixed/{{benchmark_id}}.json",
    params:
        strategy=phase3_strategy,
        pairs=int(PHASE3["fixed_effort_pairs"]),
        seed=int(PHASE3["fixed_effort_seed"]),
    log:
        f"{WORK}/logs/phase3_benchmark/inputs/{{benchmark_id}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_benchmark/inputs/{{benchmark_id}}.tsv",
    threads: 8
    resources:
        mem_mb=24000,
        runtime=180,
        partition="day",
    conda:
        "../../envs/fastp.yaml"
    script:
        "../scripts/benchmark_fixed_effort.py"


rule benchmark_reference_depleted_reads:
    input:
        r1=lambda wildcards: f"{SCRATCH}/phase3_benchmark/trimmed/{phase3_sample(wildcards)}_R1.fastq.gz",
        r2=lambda wildcards: f"{SCRATCH}/phase3_benchmark/trimmed/{phase3_sample(wildcards)}_R2.fastq.gz",
        upstream=lambda wildcards: f"{WORK}/provenance/phase3_benchmark/trim/{phase3_sample(wildcards)}.json",
        index=phase3_reference_inputs,
    output:
        r1=temp(f"{SCRATCH}/phase3_benchmark/inputs/reference/{{benchmark_id}}_R1.fastq.gz"),
        r2=temp(f"{SCRATCH}/phase3_benchmark/inputs/reference/{{benchmark_id}}_R2.fastq.gz"),
        metrics=f"{WORK}/phase3_benchmark/inputs/reference/{{benchmark_id}}.json",
        provenance=f"{WORK}/provenance/phase3_benchmark/inputs/reference/{{benchmark_id}}.json",
    params:
        sample=phase3_sample,
        strategy=phase3_strategy,
        route=lambda wildcards: SAMPLES[phase3_sample(wildcards)]["host_route"],
        prefix=PHASE3["tool_prefixes"]["assembly"],
        index_prefix=lambda wildcards: f"{WORK}/reference_indices/{SAMPLES[phase3_sample(wildcards)]['host_route']}/host",
    log:
        f"{WORK}/logs/phase3_benchmark/inputs/{{benchmark_id}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_benchmark/inputs/{{benchmark_id}}.tsv",
    threads: PHASE3["resources"]["host_depletion"]["threads"]
    resources:
        mem_mb=PHASE3["resources"]["host_depletion"]["mem_mb"],
        runtime=PHASE3["resources"]["host_depletion"]["runtime"],
        partition=PHASE3["resources"]["host_depletion"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/benchmark_reference_deplete.py"


rule benchmark_assembly:
    input:
        r1=lambda wildcards: phase3_read_input(wildcards, 1),
        r2=lambda wildcards: phase3_read_input(wildcards, 2),
        upstream=phase3_read_provenance,
        config=str(PHASE3_CONFIG_PATH),
        panel=str(PHASE3_PANEL_PATH),
    output:
        contigs=f"{WORK}/phase3_benchmark/assemblies/{{benchmark_id}}.fasta",
        metrics=f"{WORK}/phase3_benchmark/assemblies/{{benchmark_id}}.json",
        provenance=f"{WORK}/provenance/phase3_benchmark/assemblies/{{benchmark_id}}.json",
    params:
        sample=phase3_sample,
        strategy=phase3_strategy,
        prefix=PHASE3["tool_prefixes"]["assembly"],
        min_contig=int(PHASE3["min_contig_length"]),
    log:
        f"{WORK}/logs/phase3_benchmark/assemblies/{{benchmark_id}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_benchmark/assemblies/{{benchmark_id}}.tsv",
    threads: PHASE3["resources"]["assembly"]["threads"]
    resources:
        mem_mb=PHASE3["resources"]["assembly"]["mem_mb"],
        runtime=PHASE3["resources"]["assembly"]["runtime"],
        partition=PHASE3["resources"]["assembly"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/benchmark_assemble.py"


rule benchmark_markers:
    input:
        contigs=f"{WORK}/phase3_benchmark/assemblies/{{benchmark_id}}.fasta",
        upstream=f"{WORK}/provenance/phase3_benchmark/assemblies/{{benchmark_id}}.json",
    output:
        gff=f"{WORK}/phase3_benchmark/markers/{{benchmark_id}}.gff",
        metrics=f"{WORK}/phase3_benchmark/markers/{{benchmark_id}}.json",
    params:
        prefix=PHASE3["tool_prefixes"]["barrnap"],
    log:
        f"{WORK}/logs/phase3_benchmark/markers/{{benchmark_id}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_benchmark/markers/{{benchmark_id}}.tsv",
    threads: PHASE3["resources"]["markers"]["threads"]
    resources:
        mem_mb=PHASE3["resources"]["markers"]["mem_mb"],
        runtime=PHASE3["resources"]["markers"]["runtime"],
        partition=PHASE3["resources"]["markers"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/benchmark_markers.py"


rule benchmark_viruses:
    input:
        contigs=f"{WORK}/phase3_benchmark/assemblies/{{benchmark_id}}.fasta",
        upstream=f"{WORK}/provenance/phase3_benchmark/assemblies/{{benchmark_id}}.json",
        database=lambda wildcards: PHASE3["databases"]["genomad"],
    output:
        viruses=f"{WORK}/phase3_benchmark/viruses/{{benchmark_id}}.fna",
        summary=f"{WORK}/phase3_benchmark/viruses/{{benchmark_id}}.tsv",
        metrics=f"{WORK}/phase3_benchmark/viruses/{{benchmark_id}}.json",
    params:
        prefix=PHASE3["tool_prefixes"]["genomad"],
    log:
        f"{WORK}/logs/phase3_benchmark/viruses/{{benchmark_id}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_benchmark/viruses/{{benchmark_id}}.tsv",
    threads: PHASE3["resources"]["viruses"]["threads"]
    resources:
        mem_mb=PHASE3["resources"]["viruses"]["mem_mb"],
        runtime=PHASE3["resources"]["viruses"]["runtime"],
        partition=PHASE3["resources"]["viruses"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/benchmark_viruses.py"


rule benchmark_binning:
    input:
        contigs=f"{WORK}/phase3_benchmark/assemblies/{{benchmark_id}}.fasta",
        r1=lambda wildcards: phase3_read_input(wildcards, 1),
        r2=lambda wildcards: phase3_read_input(wildcards, 2),
        upstream=f"{WORK}/provenance/phase3_benchmark/assemblies/{{benchmark_id}}.json",
    output:
        bins=f"{WORK}/phase3_benchmark/bins/{{benchmark_id}}.tar.gz",
        depth=f"{WORK}/phase3_benchmark/bins/{{benchmark_id}}.depth.tsv",
        metrics=f"{WORK}/phase3_benchmark/bins/{{benchmark_id}}.json",
    params:
        prefix=PHASE3["tool_prefixes"]["assembly"],
        min_contig=int(PHASE3["metabat_min_contig"]),
        seed=int(PHASE3["metabat_seed"]),
    log:
        f"{WORK}/logs/phase3_benchmark/bins/{{benchmark_id}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_benchmark/bins/{{benchmark_id}}.tsv",
    threads: PHASE3["resources"]["binning"]["threads"]
    resources:
        mem_mb=PHASE3["resources"]["binning"]["mem_mb"],
        runtime=PHASE3["resources"]["binning"]["runtime"],
        partition=PHASE3["resources"]["binning"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/benchmark_binning.py"


rule benchmark_checkm2:
    input:
        bins=f"{WORK}/phase3_benchmark/bins/{{benchmark_id}}.tar.gz",
        upstream=f"{WORK}/phase3_benchmark/bins/{{benchmark_id}}.json",
        database=lambda wildcards: PHASE3["databases"]["checkm2"],
    output:
        quality=f"{WORK}/phase3_benchmark/checkm2/{{benchmark_id}}.tsv",
        metrics=f"{WORK}/phase3_benchmark/checkm2/{{benchmark_id}}.json",
    params:
        prefix=PHASE3["tool_prefixes"]["checkm2"],
    log:
        f"{WORK}/logs/phase3_benchmark/checkm2/{{benchmark_id}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_benchmark/checkm2/{{benchmark_id}}.tsv",
    threads: PHASE3["resources"]["checkm2"]["threads"]
    resources:
        mem_mb=PHASE3["resources"]["checkm2"]["mem_mb"],
        runtime=PHASE3["resources"]["checkm2"]["runtime"],
        partition=PHASE3["resources"]["checkm2"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/benchmark_checkm2.py"


rule benchmark_host_carryover:
    input:
        r1=lambda wildcards: phase3_read_input(wildcards, 1),
        r2=lambda wildcards: phase3_read_input(wildcards, 2),
        upstream=phase3_read_provenance,
        index=phase3_reference_inputs,
    output:
        metrics=f"{WORK}/phase3_benchmark/host_carryover/{{benchmark_id}}.json",
    params:
        sample=phase3_sample,
        route=lambda wildcards: SAMPLES[phase3_sample(wildcards)]["host_route"],
        prefix=PHASE3["tool_prefixes"]["assembly"],
        index_prefix=lambda wildcards: f"{WORK}/reference_indices/{SAMPLES[phase3_sample(wildcards)]['host_route']}/host",
    log:
        f"{WORK}/logs/phase3_benchmark/host_carryover/{{benchmark_id}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_benchmark/host_carryover/{{benchmark_id}}.tsv",
    threads: PHASE3["resources"]["host_carryover"]["threads"]
    resources:
        mem_mb=PHASE3["resources"]["host_carryover"]["mem_mb"],
        runtime=PHASE3["resources"]["host_carryover"]["runtime"],
        partition=PHASE3["resources"]["host_carryover"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/benchmark_host_carryover.py"


rule aggregate_phase3_benchmark:
    input:
        panel=str(PHASE3_PANEL_PATH),
        config=str(PHASE3_CONFIG_PATH),
        assemblies=phase3_metric_paths("assemblies"),
        markers=phase3_metric_paths("markers"),
        viruses=phase3_metric_paths("viruses"),
        bins=phase3_metric_paths("bins"),
        checkm2=phase3_metric_paths("checkm2"),
        host=phase3_metric_paths("host_carryover"),
        assembly_bench=expand(f"{WORK}/benchmarks/phase3_benchmark/assemblies/{{benchmark_id}}.tsv", benchmark_id=PHASE3_IDS),
        virus_bench=expand(f"{WORK}/benchmarks/phase3_benchmark/viruses/{{benchmark_id}}.tsv", benchmark_id=PHASE3_IDS),
        bin_bench=expand(f"{WORK}/benchmarks/phase3_benchmark/bins/{{benchmark_id}}.tsv", benchmark_id=PHASE3_IDS),
        checkm2_bench=expand(f"{WORK}/benchmarks/phase3_benchmark/checkm2/{{benchmark_id}}.tsv", benchmark_id=PHASE3_IDS),
    output:
        table=f"{WORK}/phase3_benchmark/decision_table.tsv",
        provenance=f"{WORK}/provenance/phase3_benchmark/decision_table.json",
    log:
        f"{WORK}/logs/phase3_benchmark/aggregate.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_benchmark/aggregate.tsv",
    threads: 1
    resources:
        mem_mb=PHASE3["resources"]["validation"]["mem_mb"],
        runtime=PHASE3["resources"]["validation"]["runtime"],
        partition=PHASE3["resources"]["validation"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/aggregate_phase3_benchmark.py"


rule phase3_benchmark:
    input:
        table=f"{WORK}/phase3_benchmark/decision_table.tsv",
        provenance=f"{WORK}/provenance/phase3_benchmark/decision_table.json",
    output:
        f"{WORK}/stages/phase3_benchmark.done"
    params:
        expected_rows=len(PHASE3_IDS),
        expected_samples=len(PHASE3_SAMPLES),
    log:
        f"{WORK}/logs/phase3_benchmark/validate.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_benchmark/validate.tsv",
    threads: 1
    resources:
        mem_mb=PHASE3["resources"]["validation"]["mem_mb"],
        runtime=PHASE3["resources"]["validation"]["runtime"],
        partition=PHASE3["resources"]["validation"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    script:
        "../scripts/validate_phase3_benchmark.py"
