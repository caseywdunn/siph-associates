"""
Map every library competitively against the frozen catalogs (config catalog_version) under identical
filters (config/phase5_mapping.json). Reference-bearing libraries use their
Phase-3 host-depleted reads; reference-free libraries use full capped trimmed
reads regenerated with the locked Phase-2 settings. Coverage is reported
continuously; the presence rule is locked separately from control behaviour.

1. Regenerate trimmed reads for reference-free libraries: trim_reads_for_mapping.
2. Index the combined viral catalog: index_viral_catalog_bwa.
3. Map reads and keep filtered alignments: map_reads_bacterial_bwa,
   map_reads_viral_bwa.
4. Summarize breadth and depth: summarize_bacterial_coverage_coverm,
   summarize_viral_coverage_coverm.
5. Tabulate libraries and targets, then check completeness and accounting:
   aggregate_mapping, validate_mapping.

Next: the presence rule is locked from MAG positives and decoy negatives, then
applied in Phase 6.
"""

P5 = f"{WORK}/phase5_mapping"
P5_VERSION = PHASE5["catalog_version"]
P5V = f"{P5}/{P5_VERSION}"  # catalog-dependent outputs; earlier versions stay on disk
P5_CATALOG = f"{WORK}/phase4_catalog/{PHASE5['catalog_version']}"
P5_TOOLS = f"{PHASE5['tool_prefixes']['mapping']}/bin"
P5_MAP = PHASE5["mapping"]
P5_COVERAGE = PHASE5["coverage"]
BACTERIAL_INDEX = multiext(f"{P5_CATALOG}/bacterial_catalog.fna", ".amb", ".ann", ".bwt", ".pac", ".sa")

wildcard_constraints:
    scope="pilot|cohort",


rule trim_reads_for_mapping:
    input:
        r1=raw_r1,
        r2=raw_r2,
        run_snapshot=ancient(f"{WORK}/provenance/run/inputs.sha256"),
    output:
        r1=temp(f"{P5}/trimmed/{{sample}}_R1.fastq.gz"),
        r2=temp(f"{P5}/trimmed/{{sample}}_R2.fastq.gz"),
        fastp_json=f"{P5}/trim/{{sample}}.fastp.json",
        fastp_html=f"{P5}/trim/{{sample}}.fastp.html",
        provenance=f"{WORK}/provenance/phase5_mapping/trim/{{sample}}.json",
    params:
        expected_pairs=lambda wildcards: int(SAMPLES[wildcards.sample]["read_pairs"]),
        cap_pairs=int(config["parameters"]["cap_pairs"]),
        cap_seed=int(config["parameters"]["cap_seed"]),
    log:
        f"{WORK}/logs/phase5_mapping/trim/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase5_mapping/trim/{{sample}}.tsv",
    threads: PHASE5["resources"]["trim"]["threads"]
    resources:
        **phase5_resources("trim"),
        cap_slots=lambda wildcards: 1 if int(SAMPLES[wildcards.sample]["read_pairs"]) > int(config["parameters"]["cap_pairs"]) else 0,
    conda:
        "../../envs/fastp.yaml"
    script:
        # The accepted Phase-2 trimming implementation, reused unchanged so the
        # regenerated reads match the screened reads.
        "../scripts/trim_reads.py"


rule index_viral_catalog_bwa:
    input:
        associate=f"{P5_CATALOG}/viral_associate.fna",
        endogenous=f"{P5_CATALOG}/viral_endogenous_candidate.fna",
    output:
        fasta=f"{P5V}/viral_catalog.fna",
        index=multiext(f"{P5V}/viral_catalog.fna", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    params:
        tools=P5_TOOLS,
    log:
        f"{WORK}/logs/phase5_mapping/index_viral_catalog_bwa.log",
    threads: 1
    resources:
        **phase5_resources("index"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        # Associate and endogenous-candidate vOTUs compete for the same reads.
        cat {input.associate:q} {input.endogenous:q} > {output.fasta:q}
        bwa index {output.fasta:q} > {log:q} 2>&1
        """


rule map_reads_bacterial_bwa:
    input:
        r1=lambda wildcards: phase5_reads(wildcards, 1),
        r2=lambda wildcards: phase5_reads(wildcards, 2),
        index=BACTERIAL_INDEX,
    output:
        bam=f"{P5V}/bam/bacterial/{{sample}}.bam",
        bai=f"{P5V}/bam/bacterial/{{sample}}.bam.bai",
        flagstat=f"{P5V}/bam/bacterial/{{sample}}.flagstat.txt",
    params:
        tools=P5_TOOLS,
        index=lambda wildcards, input: input.index[0].removesuffix(".amb"),
        scratch=f"{SCRATCH}/phase5_mapping/sort_tmp",
        min_score=P5_MAP["bwa_min_score"],
        min_mapq=P5_MAP["min_mapq"],
        exclude=P5_MAP["exclude_flags"],
    log:
        f"{WORK}/logs/phase5_mapping/bacterial/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase5_mapping/bacterial/{{sample}}.tsv",
    threads: PHASE5["resources"]["mapping_small"]["threads"]
    priority: 10
    resources:
        mem_mb=phase5_mapping_resource("mem_mb"),
        runtime=phase5_mapping_resource("runtime"),
        partition=phase5_mapping_resource("partition"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        mkdir -p {params.scratch:q}
        tmp=$(mktemp -d {params.scratch:q}/{wildcards.sample}.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        # Primary, mapped, properly paired alignments at MAPQ >= 30, duplicates removed.
        bwa mem -M -T {params.min_score} -t {threads} {params.index:q} {input.r1:q} {input.r2:q} 2> {log:q} \
          | samtools fixmate -m -O bam - - \
          | samtools view -u -q {params.min_mapq} -F {params.exclude} -f 0x2 - \
          | samtools sort -l 0 -@ 4 -m 1G -T "$tmp/sort" - \
          | samtools markdup -r - "$tmp/out.bam" 2>> {log:q}
        samtools index "$tmp/out.bam"
        samtools flagstat "$tmp/out.bam" > {output.flagstat:q}
        mv "$tmp/out.bam.bai" {output.bai:q}
        mv "$tmp/out.bam" {output.bam:q}
        """


rule map_reads_viral_bwa:
    input:
        r1=lambda wildcards: phase5_reads(wildcards, 1),
        r2=lambda wildcards: phase5_reads(wildcards, 2),
        index=multiext(f"{P5V}/viral_catalog.fna", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        bam=f"{P5V}/bam/viral/{{sample}}.bam",
        bai=f"{P5V}/bam/viral/{{sample}}.bam.bai",
        flagstat=f"{P5V}/bam/viral/{{sample}}.flagstat.txt",
    params:
        tools=P5_TOOLS,
        index=lambda wildcards, input: input.index[0].removesuffix(".amb"),
        scratch=f"{SCRATCH}/phase5_mapping/sort_tmp",
        min_score=P5_MAP["bwa_min_score"],
        min_mapq=P5_MAP["min_mapq"],
        exclude=P5_MAP["exclude_flags"],
    log:
        f"{WORK}/logs/phase5_mapping/viral/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase5_mapping/viral/{{sample}}.tsv",
    threads: PHASE5["resources"]["mapping_small"]["threads"]
    priority: 10
    resources:
        mem_mb=phase5_mapping_resource("mem_mb"),
        runtime=phase5_mapping_resource("runtime"),
        partition=phase5_mapping_resource("partition"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        mkdir -p {params.scratch:q}
        tmp=$(mktemp -d {params.scratch:q}/{wildcards.sample}.viral.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        bwa mem -M -T {params.min_score} -t {threads} {params.index:q} {input.r1:q} {input.r2:q} 2> {log:q} \
          | samtools fixmate -m -O bam - - \
          | samtools view -u -q {params.min_mapq} -F {params.exclude} -f 0x2 - \
          | samtools sort -l 0 -@ 4 -m 1G -T "$tmp/sort" - \
          | samtools markdup -r - "$tmp/out.bam" 2>> {log:q}
        samtools index "$tmp/out.bam"
        samtools flagstat "$tmp/out.bam" > {output.flagstat:q}
        mv "$tmp/out.bam.bai" {output.bai:q}
        mv "$tmp/out.bam" {output.bam:q}
        """


rule summarize_bacterial_coverage_coverm:
    input:
        bam=f"{P5V}/bam/bacterial/{{sample}}.bam",
        bai=f"{P5V}/bam/bacterial/{{sample}}.bam.bai",
    output:
        f"{P5V}/coverage/bacterial/{{sample}}.tsv",
    params:
        tools=P5_TOOLS,
        identity=P5_COVERAGE["min_read_percent_identity"],
        aligned=P5_COVERAGE["min_read_aligned_percent"],
        separator=P5_COVERAGE["separator"],
        methods=" ".join(P5_COVERAGE["methods"]),
    log:
        f"{WORK}/logs/phase5_mapping/coverage_bacterial/{{sample}}.log",
    threads: 4
    resources:
        **phase5_resources("index"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        coverm genome -b {input.bam:q} --separator {params.separator:q} -t {threads} \
          --min-read-percent-identity {params.identity} --min-read-aligned-percent {params.aligned} \
          --min-covered-fraction 0 -m {params.methods} > {output:q}.tmp 2> {log:q}
        mv {output:q}.tmp {output:q}
        """


rule summarize_viral_coverage_coverm:
    input:
        bam=f"{P5V}/bam/viral/{{sample}}.bam",
        bai=f"{P5V}/bam/viral/{{sample}}.bam.bai",
    output:
        f"{P5V}/coverage/viral/{{sample}}.tsv",
    params:
        tools=P5_TOOLS,
        identity=P5_COVERAGE["min_read_percent_identity"],
        aligned=P5_COVERAGE["min_read_aligned_percent"],
        methods=" ".join(P5_COVERAGE["methods"]),
    log:
        f"{WORK}/logs/phase5_mapping/coverage_viral/{{sample}}.log",
    threads: 4
    resources:
        **phase5_resources("index"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        coverm contig -b {input.bam:q} -t {threads} \
          --min-read-percent-identity {params.identity} --min-read-aligned-percent {params.aligned} \
          -m {params.methods} > {output:q}.tmp 2> {log:q}
        mv {output:q}.tmp {output:q}
        """


rule aggregate_mapping:
    input:
        script="workflow/scripts/aggregate_mapping.py",
        bacterial=lambda wildcards: expand(f"{P5V}/coverage/bacterial/{{sample}}.tsv", sample=phase5_scope_samples(wildcards)),
        viral=lambda wildcards: expand(f"{P5V}/coverage/viral/{{sample}}.tsv", sample=phase5_scope_samples(wildcards)),
        flagstat=lambda wildcards: expand(f"{P5V}/bam/{{kind}}/{{sample}}.flagstat.txt", kind=("bacterial", "viral"),
                                          sample=phase5_scope_samples(wildcards)),
        trim=lambda wildcards: [f"{WORK}/provenance/phase5_mapping/trim/{s}.json"
                                for s in phase5_scope_samples(wildcards) if SAMPLES[s]["host_route"] == "none"],
        manifest=f"{P5_CATALOG}/bacterial_catalog.manifest.tsv",
    output:
        bacterial=f"{P5V}/{{scope}}/bacterial_coverage.tsv",
        viral=f"{P5V}/{{scope}}/viral_coverage.tsv",
        libraries=f"{P5V}/{{scope}}/library_mapping.tsv",
    params:
        samples=phase5_scope_samples,
        results=lambda wildcards, output: str(Path(output.bacterial).parent.parent),
        provenance=f"{WORK}/provenance/phase5_mapping",
        eligibility=f"{WORK}/phase3_cohort/input_eligibility.tsv",
        samples_table=config["samples"],
    log:
        f"{WORK}/logs/phase5_mapping/aggregate_{{scope}}.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --samples {params.samples:q} --results {params.results:q} \
          --provenance {params.provenance:q} --manifest {input.manifest:q} --eligibility {params.eligibility:q} --samples-table {params.samples_table:q} \
          --bacterial {output.bacterial:q} --viral {output.viral:q} --libraries {output.libraries:q} > {log:q} 2>&1
        """


rule validate_mapping:
    input:
        script="workflow/scripts/validate_mapping.py",
        bacterial=f"{P5V}/{{scope}}/bacterial_coverage.tsv",
        viral=f"{P5V}/{{scope}}/viral_coverage.tsv",
        libraries=f"{P5V}/{{scope}}/library_mapping.tsv",
        manifest=f"{P5_CATALOG}/bacterial_catalog.manifest.tsv",
        viral_fasta=f"{P5V}/viral_catalog.fna",
    output:
        f"{WORK}/stages/phase5_mapping_{{scope}}_{P5_VERSION}.done",
    params:
        samples=phase5_scope_samples,
        library_qc=f"{WORK}/aggregation/cohort/library_qc.tsv",
    log:
        f"{WORK}/logs/phase5_mapping/validate_{{scope}}.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --samples {params.samples:q} --bacterial {input.bacterial:q} \
          --viral {input.viral:q} --libraries {input.libraries:q} --manifest {input.manifest:q} \
          --viral-fasta {input.viral_fasta:q} --library-qc {params.library_qc:q} --output {output:q} > {log:q} 2>&1
        """


rule phase5_pilot:
    input:
        f"{WORK}/stages/phase5_mapping_pilot_{P5_VERSION}.done",


rule phase5_mapping:
    input:
        f"{WORK}/stages/phase5_mapping_cohort_{P5_VERSION}.done",
