"""
Grade presences and run the pre-specified Phase-6 analyses
(docs/phase6_analysis_plan.md, config/phase6_analysis.json) on the frozen
Phase-5 coverage tables. Every reported table regenerates from these rules.

1. Assembly support for high-confidence grades: index_catalog_minimap2,
   align_assembly_catalog_minimap2.
2. Sensitivity read sets and mappings: trim_reads_uncapped,
   map_reads_sensitivity_bwa, summarize_sensitivity_coverage_coverm.
3. Evidence grades and physical-flowcell metadata: grade_presences,
   normalize_phase6_metadata; contamination flags: test_contamination_flowcell,
   test_nanomia_contamination_flowcell.
4. Primary analyses: summarize_incidence, fit_physalia_models_lme4 (mixed models as
   pre-specified, plus the approved within-flowcell permutation test of region).
5. Full 5% and 20% breadth analyses: analyze_presence_threshold;
   read-handling sensitivity comparisons: compare_sensitivity_analyses.
6. Check completeness of the analysis outputs: validate_phase6.

Next: Phase 7 produces publication figures and source-data tables from these outputs.
"""

P6 = f"{WORK}/phase6_analysis"
P6V = f"{P6}/{PHASE5['catalog_version']}"  # catalog-dependent outputs
P6_TOOLS = f"{PHASE5['tool_prefixes']['mapping']}/bin"
P6_STATS = f"{PHASE6['software']['stats_prefix']}/bin"
P6_CATALOG = f"{WORK}/phase4_catalog/{PHASE5['catalog_version']}"
P6_COHORT = f"{WORK}/phase5_mapping/{PHASE5['catalog_version']}/cohort"

wildcard_constraints:
    variant="full_reads|uncapped",
    kind="bacterial|viral",


rule index_catalog_minimap2:
    input:
        f"{P6_CATALOG}/bacterial_catalog.fna",
    output:
        f"{P6V}/assembly_support/bacterial_catalog.asm5.mmi",
    params:
        tools=P6_TOOLS,
    log:
        f"{WORK}/logs/phase6_analysis/index_catalog_minimap2.log",
    threads: 8
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        minimap2 -x asm5 -t {threads} -d {output:q}.tmp {input:q} > {log:q} 2>&1
        mv {output:q}.tmp {output:q}
        """


rule align_assembly_catalog_minimap2:
    input:
        contigs=f"{WORK}/phase3_cohort/assemblies/{{sample}}.fasta",
        index=f"{P6V}/assembly_support/bacterial_catalog.asm5.mmi",
    output:
        f"{P6V}/assembly_support/{{sample}}.paf",
    params:
        tools=P6_TOOLS,
    log:
        f"{WORK}/logs/phase6_analysis/assembly_support/{{sample}}.log",
    threads: 8
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        minimap2 -c --secondary=no -t {threads} {input.index:q} {input.contigs:q} > {output:q}.tmp 2> {log:q}
        mv {output:q}.tmp {output:q}
        """


rule trim_reads_uncapped:
    input:
        r1=raw_r1,
        r2=raw_r2,
        run_snapshot=ancient(f"{WORK}/provenance/run/inputs.sha256"),
    output:
        r1=temp(f"{P6}/trimmed_uncapped/{{sample}}_R1.fastq.gz"),
        r2=temp(f"{P6}/trimmed_uncapped/{{sample}}_R2.fastq.gz"),
        fastp_json=f"{P6}/trim_uncapped/{{sample}}.fastp.json",
        fastp_html=f"{P6}/trim_uncapped/{{sample}}.fastp.html",
        provenance=f"{WORK}/provenance/phase6_analysis/trim_uncapped/{{sample}}.json",
    params:
        expected_pairs=lambda wildcards: int(SAMPLES[wildcards.sample]["read_pairs"]),
        # A cap above every library's depth disables subsampling; fastp settings are unchanged.
        cap_pairs=10**12,
        cap_seed=int(config["parameters"]["cap_seed"]),
    log:
        f"{WORK}/logs/phase6_analysis/trim_uncapped/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase6_analysis/trim_uncapped/{{sample}}.tsv",
    threads: PHASE5["resources"]["trim"]["threads"]
    resources:
        **phase5_resources("trim"),
    conda:
        "../../envs/fastp.yaml"
    script:
        # The accepted Phase-2 trimming implementation, reused with capping disabled.
        "../scripts/trim_reads.py"


rule map_reads_sensitivity_bwa:
    input:
        r1=lambda wildcards: phase6_sensitivity_reads(wildcards, 1),
        r2=lambda wildcards: phase6_sensitivity_reads(wildcards, 2),
        index=phase6_index,
    output:
        bam=f"{P6V}/sensitivity/{{variant}}/bam/{{kind}}/{{sample}}.bam",
        bai=f"{P6V}/sensitivity/{{variant}}/bam/{{kind}}/{{sample}}.bam.bai",
    params:
        tools=P6_TOOLS,
        index=lambda wildcards, input: input.index[0].removesuffix(".amb"),
        scratch=f"{SCRATCH}/phase6_analysis/sort_tmp",
        min_score=PHASE5["mapping"]["bwa_min_score"],
        min_mapq=PHASE5["mapping"]["min_mapq"],
        exclude=PHASE5["mapping"]["exclude_flags"],
    log:
        f"{WORK}/logs/phase6_analysis/sensitivity/{{variant}}.{{kind}}.{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase6_analysis/sensitivity/{{variant}}.{{kind}}.{{sample}}.tsv",
    threads: 16
    priority: 10
    resources:
        **phase5_resources("mapping_large"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        mkdir -p {params.scratch:q}
        tmp=$(mktemp -d {params.scratch:q}/{wildcards.variant}.{wildcards.kind}.{wildcards.sample}.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        # Identical filters to Phase 5.
        bwa mem -M -T {params.min_score} -t {threads} {params.index:q} {input.r1:q} {input.r2:q} 2> {log:q} \
          | samtools fixmate -m -O bam - - \
          | samtools view -u -q {params.min_mapq} -F {params.exclude} -f 0x2 - \
          | samtools sort -l 0 -@ 4 -m 1G -T "$tmp/sort" - \
          | samtools markdup -r - "$tmp/out.bam" 2>> {log:q}
        samtools index "$tmp/out.bam"
        mv "$tmp/out.bam.bai" {output.bai:q}
        mv "$tmp/out.bam" {output.bam:q}
        """


rule summarize_sensitivity_coverage_coverm:
    input:
        bam=f"{P6V}/sensitivity/{{variant}}/bam/{{kind}}/{{sample}}.bam",
        bai=f"{P6V}/sensitivity/{{variant}}/bam/{{kind}}/{{sample}}.bam.bai",
    output:
        f"{P6V}/sensitivity/{{variant}}/coverage/{{kind}}/{{sample}}.tsv",
    params:
        tools=P6_TOOLS,
        mode=lambda wildcards: "genome --separator '|' --min-covered-fraction 0" if wildcards.kind == "bacterial" else "contig",
        identity=PHASE5["coverage"]["min_read_percent_identity"],
        aligned=PHASE5["coverage"]["min_read_aligned_percent"],
        methods=" ".join(PHASE5["coverage"]["methods"]),
    log:
        f"{WORK}/logs/phase6_analysis/sensitivity/{{variant}}.{{kind}}.{{sample}}.coverm.log",
    threads: 4
    resources:
        **phase5_resources("index"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        coverm {params.mode} -b {input.bam:q} -t {threads} \
          --min-read-percent-identity {params.identity} --min-read-aligned-percent {params.aligned} \
          -m {params.methods} > {output:q}.tmp 2> {log:q}
        mv {output:q}.tmp {output:q}
        """


rule grade_presences:
    input:
        script="workflow/scripts/grade_presences.py",
        presence=str(ROOT / "config" / "phase5_presence.json"),
        analysis=str(PHASE6_CONFIG_PATH),
        bacterial=f"{P6_COHORT}/bacterial_coverage.tsv",
        viral=f"{P6_COHORT}/viral_coverage.tsv",
        manifest=f"{P6_CATALOG}/bacterial_catalog.manifest.tsv",
        nominated=f"{WORK}/phase4_catalog/references/nominated.tsv",
        mags=f"{WORK}/phase3_catalog/mags/mag_catalog.tsv",
        votus=[f"{WORK}/phase3_catalog/votus/{c}.votus.tsv" for c in ("associate", "endogenous_candidate")],
        support=expand(f"{P6V}/assembly_support/{{sample}}.paf", sample=ASSEMBLY_IDS),
    output:
        bacterial=f"{P6V}/grades/bacterial_grades.tsv",
        viral=f"{P6V}/grades/viral_grades.tsv",
    params:
        support_dir=lambda wildcards, input: str(Path(input.support[0]).parent),
    log:
        f"{WORK}/logs/phase6_analysis/grade_presences.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --presence {input.presence:q} --analysis {input.analysis:q} \
          --bacterial {input.bacterial:q} --viral {input.viral:q} --manifest {input.manifest:q} \
          --nominated {input.nominated:q} --mags {input.mags:q} --votus {input.votus:q} \
          --support-dir {params.support_dir:q} --out-bacterial {output.bacterial:q} \
          --out-viral {output.viral:q} > {log:q} 2>&1
        """


rule normalize_phase6_metadata:
    input:
        script="workflow/scripts/phase6_metadata.py",
        manifest=config["manifest"],
    output:
        f"{P6V}/primary/library_flowcells.tsv",
    log:
        f"{WORK}/logs/phase6_analysis/normalize_phase6_metadata.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --manifest {input.manifest:q} --flowcell-policy single_flowcell \
          --output {output:q} > {log:q} 2>&1
        """


rule test_contamination_flowcell:
    input:
        script="workflow/scripts/test_contamination_flowcell.py",
        analysis=str(PHASE6_CONFIG_PATH),
        grades=f"{P6V}/grades/bacterial_grades.tsv",
        metadata=f"{P6V}/primary/library_flowcells.tsv",
    output:
        f"{P6V}/grades/contamination_tests.tsv",
    log:
        f"{WORK}/logs/phase6_analysis/test_contamination_flowcell.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --analysis {input.analysis:q} --grades {input.grades:q} \
          --metadata {input.metadata:q} --output {output:q} > {log:q} 2>&1
        """


rule test_nanomia_contamination_flowcell:
    input:
        script="workflow/scripts/test_contamination_flowcell.py",
        analysis=str(PHASE6_CONFIG_PATH),
        grades=f"{P6V}/grades/bacterial_grades.tsv",
        metadata=f"{P6V}/primary/library_flowcells.tsv",
    output:
        f"{P6V}/primary/nanomia_contamination_tests.tsv",
    log:
        f"{WORK}/logs/phase6_analysis/test_nanomia_contamination_flowcell.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --analysis {input.analysis:q} --grades {input.grades:q} \
          --metadata {input.metadata:q} --subset nanomia --output {output:q} > {log:q} 2>&1
        """


rule summarize_incidence:
    input:
        script="workflow/scripts/summarize_incidence.py",
        bacterial=f"{P6V}/grades/bacterial_grades.tsv",
        viral=f"{P6V}/grades/viral_grades.tsv",
        contamination=f"{P6V}/grades/contamination_tests.tsv",
        manifest=config["manifest"],
        samples=config["samples"],
        links=f"{P6_CATALOG}/crispr/host_links.tsv",
    output:
        bacterial=f"{P6V}/primary/bacterial_incidence.tsv",
        viral=f"{P6V}/primary/viral_incidence.tsv",
        grades=f"{P6V}/primary/grade_summary.tsv",
        links=f"{P6V}/primary/phage_host_links.tsv",
    log:
        f"{WORK}/logs/phase6_analysis/summarize_incidence.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --bacterial {input.bacterial:q} --viral {input.viral:q} \
          --contamination {input.contamination:q} --manifest {input.manifest:q} --samples {input.samples:q} \
          --links {input.links:q} \
          --out-bacterial {output.bacterial:q} --out-viral {output.viral:q} --out-grades {output.grades:q} \
          --out-links {output.links:q} > {log:q} 2>&1
        """


rule fit_physalia_models_lme4:
    input:
        script="workflow/scripts/fit_physalia_models.R",
        analysis=str(PHASE6_CONFIG_PATH),
        grades=f"{P6V}/grades/bacterial_grades.tsv",
        libraries=f"{P6_COHORT}/library_mapping.tsv",
        metadata=f"{P6V}/primary/library_flowcells.tsv",
    output:
        design=f"{P6V}/primary/physalia_design.tsv",
        models=f"{P6V}/primary/physalia_models.tsv",
        permanova=f"{P6V}/primary/physalia_permanova.tsv",
        permutation=f"{P6V}/primary/physalia_region_permutation.tsv",
    params:
        stats=P6_STATS,
    log:
        f"{WORK}/logs/phase6_analysis/fit_physalia_models.log",
    threads: 4
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        {params.stats}/Rscript {input.script:q} {input.analysis:q} {input.grades:q} {input.libraries:q} \
          {input.metadata:q} {output.design:q} {output.models:q} {output.permanova:q} {output.permutation:q} grade {threads} > {log:q} 2>&1
        """


rule compare_sensitivity_analyses:
    input:
        script="workflow/scripts/compare_sensitivity.py",
        presence=str(ROOT / "config" / "phase5_presence.json"),
        subsets=str(ROOT / "config" / "phase6_subsets.tsv"),
        grades=f"{P6V}/grades/bacterial_grades.tsv",
        cohort_viral=f"{P6_COHORT}/viral_coverage.tsv",
        nominated=f"{WORK}/phase4_catalog/references/nominated.tsv",
        catalog=f"{P6_CATALOG}/bacterial_catalog.manifest.tsv",
        manifest=config["manifest"],
        coverage=[f"{P6V}/sensitivity/{v}/coverage/{k}/{s}.tsv" for v, s in PHASE6_SENSITIVITY
                  for k in ("bacterial", "viral")],
    output:
        thresholds=f"{P6V}/sensitivity/threshold_sensitivity.tsv",
        host_handling=f"{P6V}/sensitivity/host_handling_concordance.tsv",
        capping=f"{P6V}/sensitivity/capping_concordance.tsv",
        loso=f"{P6V}/sensitivity/leave_one_study_out.tsv",
    params:
        sensitivity_dir=lambda wildcards, output: str(Path(output.thresholds).parent),
    log:
        f"{WORK}/logs/phase6_analysis/compare_sensitivity.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --presence {input.presence:q} --subsets {input.subsets:q} \
          --grades {input.grades:q} --cohort-viral {input.cohort_viral:q} --nominated {input.nominated:q} \
          --catalog {input.catalog:q} --manifest {input.manifest:q} --sensitivity-dir {params.sensitivity_dir:q} \
          --out-thresholds {output.thresholds:q} --out-host-handling {output.host_handling:q} \
          --out-capping {output.capping:q} --out-loso {output.loso:q} > {log:q} 2>&1
        """


rule analyze_presence_threshold:
    input:
        contamination_script="workflow/scripts/test_contamination_flowcell.py",
        incidence_script="workflow/scripts/summarize_incidence.py",
        models_script="workflow/scripts/fit_physalia_models.R",
        analysis=str(PHASE6_CONFIG_PATH),
        bacterial=f"{P6V}/grades/bacterial_grades.tsv",
        viral=f"{P6V}/grades/viral_grades.tsv",
        libraries=f"{P6_COHORT}/library_mapping.tsv",
        metadata=f"{P6V}/primary/library_flowcells.tsv",
        manifest=config["manifest"],
        samples=config["samples"],
        links=f"{P6_CATALOG}/crispr/host_links.tsv",
    output:
        contamination=f"{P6V}/sensitivity/breadth_{{threshold}}pct/contamination_tests.tsv",
        nanomia=f"{P6V}/sensitivity/breadth_{{threshold}}pct/nanomia_contamination_tests.tsv",
        bacterial=f"{P6V}/sensitivity/breadth_{{threshold}}pct/bacterial_incidence.tsv",
        viral=f"{P6V}/sensitivity/breadth_{{threshold}}pct/viral_incidence.tsv",
        grades=f"{P6V}/sensitivity/breadth_{{threshold}}pct/grade_summary.tsv",
        links=f"{P6V}/sensitivity/breadth_{{threshold}}pct/phage_host_links.tsv",
        design=f"{P6V}/sensitivity/breadth_{{threshold}}pct/physalia_design.tsv",
        models=f"{P6V}/sensitivity/breadth_{{threshold}}pct/physalia_models.tsv",
        permanova=f"{P6V}/sensitivity/breadth_{{threshold}}pct/physalia_permanova.tsv",
        permutation=f"{P6V}/sensitivity/breadth_{{threshold}}pct/physalia_region_permutation.tsv",
    params:
        stats=P6_STATS,
        grade=lambda wildcards: f"grade_at_{wildcards.threshold}pct",
    wildcard_constraints:
        threshold="5|20",
    log:
        f"{WORK}/logs/phase6_analysis/breadth_{{threshold}}pct.log",
    threads: 4
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.contamination_script:q} --analysis {input.analysis:q} --grades {input.bacterial:q} \
          --metadata {input.metadata:q} --grade-column {params.grade:q} \
          --output {output.contamination:q} > {log:q} 2>&1
        python {input.contamination_script:q} --analysis {input.analysis:q} --grades {input.bacterial:q} \
          --metadata {input.metadata:q} --grade-column {params.grade:q} --subset nanomia \
          --output {output.nanomia:q} >> {log:q} 2>&1
        python {input.incidence_script:q} --bacterial {input.bacterial:q} --viral {input.viral:q} \
          --contamination {output.contamination:q} --manifest {input.manifest:q} --samples {input.samples:q} \
          --links {input.links:q} --grade-column {params.grade:q} \
          --out-bacterial {output.bacterial:q} --out-viral {output.viral:q} --out-grades {output.grades:q} \
          --out-links {output.links:q} >> {log:q} 2>&1
        {params.stats}/Rscript {input.models_script:q} {input.analysis:q} {input.bacterial:q} {input.libraries:q} \
          {input.metadata:q} {output.design:q} {output.models:q} {output.permanova:q} \
          {output.permutation:q} {params.grade:q} {threads} >> {log:q} 2>&1
        """


rule validate_phase6:
    input:
        script="workflow/scripts/validate_phase6.py",
        bacterial=f"{P6V}/grades/bacterial_grades.tsv",
        viral=f"{P6V}/grades/viral_grades.tsv",
        contamination=f"{P6V}/grades/contamination_tests.tsv",
        metadata=f"{P6V}/primary/library_flowcells.tsv",
        nanomia=f"{P6V}/primary/nanomia_contamination_tests.tsv",
        thresholds=expand(f"{P6V}/sensitivity/breadth_{{threshold}}pct/{{name}}.tsv",
                          threshold=[5, 20], name=["contamination_tests", "nanomia_contamination_tests",
                          "bacterial_incidence", "viral_incidence", "grade_summary", "phage_host_links",
                          "physalia_design", "physalia_models", "physalia_permanova", "physalia_region_permutation"]),
        primary=[f"{P6V}/primary/{name}.tsv" for name in ("bacterial_incidence", "viral_incidence", "grade_summary",
                                                         "phage_host_links", "physalia_design", "physalia_models",
                                                         "physalia_permanova", "physalia_region_permutation")],
        sensitivity=[f"{P6V}/sensitivity/{name}.tsv" for name in ("threshold_sensitivity", "host_handling_concordance",
                                                                 "capping_concordance", "leave_one_study_out")],
    output:
        f"{WORK}/stages/phase6_analysis_{PHASE5['catalog_version']}.done",
    log:
        f"{WORK}/logs/phase6_analysis/validate_phase6.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --bacterial {input.bacterial:q} --viral {input.viral:q} \
          --contamination {input.contamination:q} --metadata {input.metadata:q} --nanomia {input.nanomia:q} \
          --thresholds {input.thresholds:q} --primary {input.primary:q} \
          --sensitivity {input.sensitivity:q} --output {output:q} > {log:q} 2>&1
        """


rule phase6_analysis:
    input:
        f"{WORK}/stages/phase6_analysis_{PHASE5['catalog_version']}.done",
