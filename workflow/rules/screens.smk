rule kraken_bracken:
    input:
        r1=ancient(f"{SCRATCH}/trimmed/{{sample}}_R1.fastq.gz"),
        r2=ancient(f"{SCRATCH}/trimmed/{{sample}}_R2.fastq.gz"),
        fastp=f"{WORK}/trim/{{sample}}.fastp.json",
        trim_provenance=f"{WORK}/provenance/trim/{{sample}}.json",
        database=lambda wildcards: config["databases"]["kraken2"],
    output:
        report=f"{WORK}/screens/kraken/{{sample}}.report.tsv",
        assignments=f"{WORK}/screens/kraken/{{sample}}.assignments.tsv.gz",
        classified_r1=f"{WORK}/screens/kraken/{{sample}}.classified_1.fastq.gz",
        classified_r2=f"{WORK}/screens/kraken/{{sample}}.classified_2.fastq.gz",
        genus=f"{WORK}/screens/bracken/{{sample}}.G.tsv",
        species=f"{WORK}/screens/bracken/{{sample}}.S.tsv",
        provenance=f"{WORK}/provenance/kraken_bracken/{{sample}}.json",
    params:
        prefix=config["tool_prefixes"]["bracken"],
        threshold=int(config["parameters"]["bracken"]["threshold"]),
        read_length=int(config["parameters"]["bracken"]["read_length"]),
    log:
        f"{WORK}/logs/kraken_bracken/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/kraken_bracken/{{sample}}.tsv",
    threads: config["resources"]["kraken_bracken"]["threads"]
    priority: 100
    resources:
        mem_mb=config["resources"]["kraken_bracken"]["mem_mb"],
        runtime=config["resources"]["kraken_bracken"]["runtime"],
        partition=config["resources"]["kraken_bracken"]["partition"],
    conda:
        "../../envs/kraken-bracken.yaml"
    script:
        "../scripts/kraken_bracken.py"


rule sylph_screen:
    input:
        r1=ancient(f"{SCRATCH}/trimmed/{{sample}}_R1.fastq.gz"),
        r2=ancient(f"{SCRATCH}/trimmed/{{sample}}_R2.fastq.gz"),
        trim_provenance=f"{WORK}/provenance/trim/{{sample}}.json",
        databases=lambda wildcards: config["databases"]["sylph"],
    output:
        profile=f"{WORK}/screens/sylph/{{sample}}.profile.tsv",
        provenance=f"{WORK}/provenance/sylph/{{sample}}.json",
    params:
        prefix=config["tool_prefixes"]["sylph"],
    log:
        f"{WORK}/logs/sylph/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/sylph/{{sample}}.tsv",
    threads: config["resources"]["sylph"]["threads"]
    priority: 100
    resources:
        mem_mb=config["resources"]["sylph"]["mem_mb"],
        runtime=config["resources"]["sylph"]["runtime"],
        partition=config["resources"]["sylph"]["partition"],
    conda:
        "../../envs/sylph.yaml"
    script:
        "../scripts/sylph_screen.py"


rule phyloflash_screen:
    input:
        r1=ancient(f"{SCRATCH}/trimmed/{{sample}}_R1.fastq.gz"),
        r2=ancient(f"{SCRATCH}/trimmed/{{sample}}_R2.fastq.gz"),
        trim_provenance=f"{WORK}/provenance/trim/{{sample}}.json",
        database=lambda wildcards: config["databases"]["phyloflash"],
    output:
        archive=f"{WORK}/screens/phyloflash/{{sample}}.tar.gz",
        provenance=f"{WORK}/provenance/phyloflash/{{sample}}.json",
    params:
        prefix=config["tool_prefixes"]["phyloflash"],
    log:
        f"{WORK}/logs/phyloflash/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/phyloflash/{{sample}}.tsv",
    threads: config["resources"]["phyloflash"]["threads"]
    priority: 100
    resources:
        mem_mb=config["resources"]["phyloflash"]["mem_mb"],
        runtime=config["resources"]["phyloflash"]["runtime"],
        partition=config["resources"]["phyloflash"]["partition"],
    # phyloFlash 3.4.2 requires its Python-2-era EMIRGE environment. Keep the
    # Python 3 workflow wrapper outside that environment; it prepends the
    # versioned phyloFlash tool prefix to the child process PATH.
    conda:
        "../../envs/phyloflash.yaml"
    script:
        "../scripts/phyloflash_screen.py"
