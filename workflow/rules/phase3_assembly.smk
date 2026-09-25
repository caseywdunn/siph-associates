"""
Assemble every eligible library and recover its markers, viruses, and MAG bins.
Inputs are the accepted Phase-3 cohort reads (reference-depleted for P_physalis
and N_septata, a fixed 25 M-pair subsample for reference-free hosts), limited to
the frozen membership in config/phase3_assembly.tsv. Products seed the pooled
catalogs; they are not abundance measures.

1. Assemble contigs per library: assemble_contigs_megahit.
2. Find rRNA markers on the contigs: find_rrna_barrnap.
3. Identify and assess viral contigs: identify_viruses_genomad, assess_viruses_checkv.
4. Bin contigs by back-mapped coverage and assess bins: bin_contigs_metabat2,
   assess_bins_checkm2.
5. Tabulate all libraries, including documented exclusions, and check exact
   membership: summarize_assemblies, validate_assemblies.

Next: pooled GTDB-Tk classification, locked MAG/vOTU quality rules, and cohort
dereplication consume these per-library products in a later stage.
"""

ASSEMBLY_RESOURCES = PHASE3_ASSEMBLY["resources"]
ASSEMBLY_TOOLS = {name: f"{prefix}/bin" for name, prefix in PHASE3_ASSEMBLY["tool_prefixes"].items()}
ASSEMBLY_SCRATCH = f"{SCRATCH}/phase3_cohort"


rule assemble_contigs_megahit:
    input:
        r1=lambda wildcards: phase3_cohort_input_path(wildcards, 1),
        r2=lambda wildcards: phase3_cohort_input_path(wildcards, 2),
    output:
        contigs=f"{WORK}/phase3_cohort/assemblies/{{sample}}.fasta",
    params:
        tools=ASSEMBLY_TOOLS["assembly"],
        scratch=f"{ASSEMBLY_SCRATCH}/megahit_tmp",
        min_contig=int(PHASE3_ASSEMBLY["min_contig_length"]),
    log:
        f"{WORK}/logs/phase3_cohort/assemblies/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_cohort/assemblies/{{sample}}.tsv",
    threads: ASSEMBLY_RESOURCES["assembly_large"]["threads"]
    resources:
        mem_mb=phase3_assembly_resource("mem_mb"),
        runtime=phase3_assembly_resource("runtime"),
        partition=phase3_assembly_resource("partition"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        mkdir -p {params.scratch:q}
        tmp=$(mktemp -d {params.scratch:q}/{wildcards.sample}.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        megahit --version > {log:q} 2>&1
        megahit -1 {input.r1:q} -2 {input.r2:q} -t {threads} -o "$tmp/assembly" \
          --min-contig-len {params.min_contig} >> {log:q} 2>&1
        # Prefix contig IDs with the library so every downstream sequence keeps its source.
        sed 's/^>/>{wildcards.sample}__/' "$tmp/assembly/final.contigs.fa" > "$tmp/contigs.fasta"
        mv "$tmp/contigs.fasta" {output.contigs:q}
        """


rule find_rrna_barrnap:
    input:
        contigs=f"{WORK}/phase3_cohort/assemblies/{{sample}}.fasta",
    output:
        gff=f"{WORK}/phase3_cohort/markers/{{sample}}.gff",
    params:
        tools=ASSEMBLY_TOOLS["barrnap"],
    log:
        f"{WORK}/logs/phase3_cohort/markers/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_cohort/markers/{{sample}}.tsv",
    threads: ASSEMBLY_RESOURCES["markers"]["threads"]
    resources:
        mem_mb=ASSEMBLY_RESOURCES["markers"]["mem_mb"],
        runtime=ASSEMBLY_RESOURCES["markers"]["runtime"],
        partition=ASSEMBLY_RESOURCES["markers"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        barrnap --version > {log:q} 2>&1
        : > {output.gff:q}.tmp
        if [ -s {input.contigs:q} ]; then
          for kingdom in bac arc euk; do
            barrnap --kingdom "$kingdom" --threads {threads} {input.contigs:q} 2>> {log:q} \
              | awk -F '\t' -v OFS='\t' -v kingdom="$kingdom" \
                  '!/^#/ && NF {{ sub(/;$/, "", $NF); $NF = $NF ";kingdom=" kingdom; print }}' \
              >> {output.gff:q}.tmp
          done
        else
          echo "Assembly has no contigs; barrnap not run." >> {log:q}
        fi
        mv {output.gff:q}.tmp {output.gff:q}
        """


rule identify_viruses_genomad:
    input:
        contigs=f"{WORK}/phase3_cohort/assemblies/{{sample}}.fasta",
        database=PHASE3_ASSEMBLY["databases"]["genomad"],
    output:
        viruses=f"{WORK}/phase3_cohort/viruses/{{sample}}.fna",
        summary=f"{WORK}/phase3_cohort/viruses/{{sample}}.tsv",
    params:
        tools=ASSEMBLY_TOOLS["genomad"],
        scratch=f"{ASSEMBLY_SCRATCH}/genomad_tmp",
    log:
        f"{WORK}/logs/phase3_cohort/viruses/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_cohort/viruses/{{sample}}.tsv",
    threads: ASSEMBLY_RESOURCES["viruses"]["threads"]
    resources:
        mem_mb=ASSEMBLY_RESOURCES["viruses"]["mem_mb"],
        runtime=ASSEMBLY_RESOURCES["viruses"]["runtime"],
        partition=ASSEMBLY_RESOURCES["viruses"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        mkdir -p {params.scratch:q}
        tmp=$(mktemp -d {params.scratch:q}/{wildcards.sample}.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        genomad --version > {log:q} 2>&1
        summary="$tmp/genomad/assembly_summary"
        if [ -s {input.contigs:q} ]; then
          ln -s "$(realpath {input.contigs:q})" "$tmp/assembly.fasta"
          genomad end-to-end --threads {threads} --cleanup \
            "$tmp/assembly.fasta" "$tmp/genomad" {input.database:q} >> {log:q} 2>&1
        else
          echo "Assembly has no contigs; geNomad not run." >> {log:q}
        fi
        mkdir -p "$summary"
        [ -f "$summary/assembly_virus.fna" ] || : > "$summary/assembly_virus.fna"
        [ -f "$summary/assembly_virus_summary.tsv" ] || echo seq_name > "$summary/assembly_virus_summary.tsv"
        mv "$summary/assembly_virus_summary.tsv" {output.summary:q}
        mv "$summary/assembly_virus.fna" {output.viruses:q}
        """


rule assess_viruses_checkv:
    input:
        viruses=f"{WORK}/phase3_cohort/viruses/{{sample}}.fna",
        database=PHASE3_ASSEMBLY["databases"]["checkv"],
    output:
        quality=f"{WORK}/phase3_cohort/checkv/{{sample}}.tsv",
    params:
        tools=ASSEMBLY_TOOLS["genomad"],
        scratch=f"{ASSEMBLY_SCRATCH}/checkv_tmp",
    log:
        f"{WORK}/logs/phase3_cohort/checkv/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_cohort/checkv/{{sample}}.tsv",
    threads: ASSEMBLY_RESOURCES["checkv"]["threads"]
    resources:
        mem_mb=ASSEMBLY_RESOURCES["checkv"]["mem_mb"],
        runtime=ASSEMBLY_RESOURCES["checkv"]["runtime"],
        partition=ASSEMBLY_RESOURCES["checkv"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        mkdir -p {params.scratch:q}
        tmp=$(mktemp -d {params.scratch:q}/{wildcards.sample}.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        checkv -h > "$tmp/help" 2>&1
        head -n 1 "$tmp/help" > {log:q}
        if [ -s {input.viruses:q} ]; then
          checkv end_to_end {input.viruses:q} "$tmp/checkv" -d {input.database:q} \
            -t {threads} --remove_tmp >> {log:q} 2>&1
          mv "$tmp/checkv/quality_summary.tsv" {output.quality:q}
        else
          echo "No geNomad viral contigs; CheckV not run." >> {log:q}
          printf 'contig_id\tcontig_length\tcheckv_quality\tcompleteness\tcontamination\n' > {output.quality:q}
        fi
        """


rule bin_contigs_metabat2:
    input:
        contigs=f"{WORK}/phase3_cohort/assemblies/{{sample}}.fasta",
        r1=lambda wildcards: phase3_cohort_input_path(wildcards, 1),
        r2=lambda wildcards: phase3_cohort_input_path(wildcards, 2),
    output:
        bins=f"{WORK}/phase3_cohort/bins/{{sample}}.tar.gz",
        depth=f"{WORK}/phase3_cohort/bins/{{sample}}.depth.tsv",
    params:
        tools=ASSEMBLY_TOOLS["assembly"],
        scratch=f"{ASSEMBLY_SCRATCH}/binning_tmp",
        min_contig=int(PHASE3_ASSEMBLY["metabat_min_contig"]),
        seed=int(PHASE3_ASSEMBLY["metabat_seed"]),
    log:
        f"{WORK}/logs/phase3_cohort/bins/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_cohort/bins/{{sample}}.tsv",
    threads: ASSEMBLY_RESOURCES["binning"]["threads"]
    resources:
        mem_mb=ASSEMBLY_RESOURCES["binning"]["mem_mb"],
        runtime=ASSEMBLY_RESOURCES["binning"]["runtime"],
        partition=ASSEMBLY_RESOURCES["binning"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        mkdir -p {params.scratch:q}
        tmp=$(mktemp -d {params.scratch:q}/{wildcards.sample}.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        mkdir "$tmp/bins"
        metabat2 --help > "$tmp/help" 2>&1
        head -n 2 "$tmp/help" > {log:q}
        binnable=$(awk -v min={params.min_contig} \
          '/^>/ {{ if (length_ >= min) n++; length_ = 0; next }} {{ length_ += length($0) }}
           END {{ if (length_ >= min) n++; print n + 0 }}' {input.contigs:q})
        if [ "$binnable" -gt 0 ]; then
          # The sorted BAM stays in scratch; coverage persists as the depth table.
          ln -s "$(realpath {input.contigs:q})" "$tmp/assembly.fasta"
          bwa index "$tmp/assembly.fasta" >> {log:q} 2>&1
          bwa mem -t {threads} "$tmp/assembly.fasta" {input.r1:q} {input.r2:q} 2>> {log:q} \
            | samtools sort -@ $(( {threads} / 2 )) -o "$tmp/reads.sorted.bam" - 2>> {log:q}
          jgi_summarize_bam_contig_depths --outputDepth "$tmp/depth.tsv" "$tmp/reads.sorted.bam" >> {log:q} 2>&1
          metabat2 -i "$tmp/assembly.fasta" -a "$tmp/depth.tsv" -o "$tmp/bins/{wildcards.sample}__bin" \
            -t {threads} -m {params.min_contig} --seed {params.seed} >> {log:q} 2>&1
        else
          echo "No contigs >= {params.min_contig} bp; back-mapping and MetaBAT2 not run." >> {log:q}
          printf 'contigName\tcontigLen\ttotalAvgDepth\n' > "$tmp/depth.tsv"
        fi
        (cd "$tmp/bins" && shopt -s nullglob && bins=(*.fa) && \
          tar -czf "$tmp/bins.tar.gz" --files-from /dev/null "${{bins[@]}}")
        mv "$tmp/depth.tsv" {output.depth:q}
        mv "$tmp/bins.tar.gz" {output.bins:q}
        """


rule assess_bins_checkm2:
    input:
        bins=f"{WORK}/phase3_cohort/bins/{{sample}}.tar.gz",
        database=PHASE3_ASSEMBLY["databases"]["checkm2"],
    output:
        quality=f"{WORK}/phase3_cohort/checkm2/{{sample}}.tsv",
        status=f"{WORK}/phase3_cohort/checkm2/{{sample}}.status",
    params:
        tools=ASSEMBLY_TOOLS["checkm2"],
        scratch=f"{ASSEMBLY_SCRATCH}/checkm2_tmp",
        no_annotations="No DIAMOND annotation was generated",
    log:
        f"{WORK}/logs/phase3_cohort/checkm2/{{sample}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_cohort/checkm2/{{sample}}.tsv",
    threads: ASSEMBLY_RESOURCES["checkm2"]["threads"]
    resources:
        mem_mb=ASSEMBLY_RESOURCES["checkm2"]["mem_mb"],
        runtime=ASSEMBLY_RESOURCES["checkm2"]["runtime"],
        partition=ASSEMBLY_RESOURCES["checkm2"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        mkdir -p {params.scratch:q}
        tmp=$(mktemp -d {params.scratch:q}/{wildcards.sample}.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        mkdir "$tmp/bins"
        tar -xzf {input.bins:q} -C "$tmp/bins"
        checkm2 --version > {log:q} 2>&1
        empty_report='Name\tCompleteness\tContamination\n'
        if ! compgen -G "$tmp/bins/*.fa" > /dev/null; then
          echo "No MetaBAT2 bins; CheckM2 not run." >> {log:q}
          status=no_bins
          printf "$empty_report" > "$tmp/quality.tsv"
        elif checkm2 predict --threads {threads} --input "$tmp/bins" --extension fa \
            --database_path {input.database:q} --output-directory "$tmp/output" --force >> {log:q} 2>&1; then
          status=assessed
          mv "$tmp/output/quality_report.tsv" "$tmp/quality.tsv"
        elif grep -qF {params.no_annotations:q} {log:q}; then
          # A distinct state from low quality: CheckM2 found nothing to assess.
          echo "Recorded bins as unassessable because CheckM2 found no DIAMOND annotations." >> {log:q}
          status=no_annotations
          printf "$empty_report" > "$tmp/quality.tsv"
        else
          exit 1
        fi
        mv "$tmp/quality.tsv" {output.quality:q}
        echo "$status" > {output.status:q}
        """


rule summarize_assemblies:
    input:
        script="workflow/scripts/summarize_assemblies.py",
        membership=str(PHASE3_ASSEMBLY_MEMBERSHIP_PATH),
        contigs=expand(f"{WORK}/phase3_cohort/assemblies/{{sample}}.fasta", sample=ASSEMBLY_IDS),
        markers=expand(f"{WORK}/phase3_cohort/markers/{{sample}}.gff", sample=ASSEMBLY_IDS),
        viruses=expand(f"{WORK}/phase3_cohort/viruses/{{sample}}.fna", sample=ASSEMBLY_IDS),
        checkv=expand(f"{WORK}/phase3_cohort/checkv/{{sample}}.tsv", sample=ASSEMBLY_IDS),
        bins=expand(f"{WORK}/phase3_cohort/bins/{{sample}}.tar.gz", sample=ASSEMBLY_IDS),
        checkm2=expand(f"{WORK}/phase3_cohort/checkm2/{{sample}}.tsv", sample=ASSEMBLY_IDS),
        checkm2_status=expand(f"{WORK}/phase3_cohort/checkm2/{{sample}}.status", sample=ASSEMBLY_IDS),
        assembly_bench=expand(f"{WORK}/benchmarks/phase3_cohort/assemblies/{{sample}}.tsv", sample=ASSEMBLY_IDS),
    output:
        table=f"{WORK}/phase3_cohort/assembly_summary.tsv",
    params:
        results=lambda wildcards, output: str(Path(output.table).parent),
        benchmarks=lambda wildcards, input: str(Path(input.assembly_bench[0]).parent),
        metabat_min_contig=int(PHASE3_ASSEMBLY["metabat_min_contig"]),
    log:
        f"{WORK}/logs/phase3_cohort/assembly_summary.log",
    threads: 1
    resources:
        mem_mb=ASSEMBLY_RESOURCES["validation"]["mem_mb"],
        runtime=ASSEMBLY_RESOURCES["validation"]["runtime"],
        partition=ASSEMBLY_RESOURCES["validation"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --membership {input.membership:q} --results {params.results:q} \
          --benchmarks {params.benchmarks:q} --metabat-min-contig {params.metabat_min_contig} \
          --output {output.table:q} > {log:q} 2>&1
        """


rule validate_assemblies:
    input:
        script="workflow/scripts/validate_assemblies.py",
        table=f"{WORK}/phase3_cohort/assembly_summary.tsv",
        membership=str(PHASE3_ASSEMBLY_MEMBERSHIP_PATH),
        eligibility=f"{WORK}/phase3_cohort/input_eligibility.tsv",
    output:
        f"{WORK}/stages/phase3_assembly.done"
    log:
        f"{WORK}/logs/phase3_cohort/validate_assembly.log",
    threads: 1
    resources:
        mem_mb=ASSEMBLY_RESOURCES["validation"]["mem_mb"],
        runtime=ASSEMBLY_RESOURCES["validation"]["runtime"],
        partition=ASSEMBLY_RESOURCES["validation"]["partition"],
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --summary {input.table:q} --membership {input.membership:q} \
          --eligibility {input.eligibility:q} --output {output:q} > {log:q} 2>&1
        """


rule phase3_assembly:
    input:
        f"{WORK}/stages/phase3_assembly.done",
