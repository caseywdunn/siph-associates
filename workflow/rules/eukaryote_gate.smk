"""
Eukaryotic parasite and prey evidence from SSU rRNA. The accepted Phase-2
phyloFlash archives supply each library's SSU reads and assembled full-length
SSUs. Reads are re-mapped competitively to SILVA 138.1 plus this cohort's host
siphonophore SSUs and assigned to the lowest common ancestor of their near-tied
hits, so host reads and reads from conserved 18S regions do not masquerade as
other animals.

1. Extract SSU reads and assembled SSUs per library: extract_phyloflash_ssu.
2. Build the competitive reference with host SSUs: build_eukaryote_ssu_reference,
   index_eukaryote_ssu_reference_minimap2.
3. Assign read pairs to lineages: classify_ssu_reads_minimap2.
4. Tabulate all libraries: aggregate_ssu_read_lineages.
5. Classify assembled SSUs from phyloFlash and from Phase-3 assemblies:
   classify_assembled_ssu_minimap2, collect_eukaryote_assembled.
6. Verify assembled non-host SSUs at NCBI: verify_eukaryotes_ncbi.
7. Grade detections under the locked rules and check them: grade_eukaryotes,
   validate_eukaryote_gate.

Rules are in config/eukaryote_gate.json; reasoning in docs/eukaryote_gate_decisions.md.
"""

EUK = f"{WORK}/eukaryote_gate"
EUK_TOOLS = f"{PHASE5['tool_prefixes']['mapping']}/bin"
EUK_SILVA = EUKARYOTE["databases"]["silva_ssu_nr99"]



rule extract_phyloflash_ssu:
    input:
        script="workflow/scripts/extract_phyloflash_ssu.py",
    output:
        r1=f"{EUK}/ssu/{{sample}}.R1.fq.gz",
        r2=f"{EUK}/ssu/{{sample}}.R2.fq.gz",
        assembled=f"{EUK}/ssu/{{sample}}.assembled_ssu.fasta",
        classification=f"{EUK}/ssu/{{sample}}.ssu_classification.csv",
    params:
        # Accepted Phase-2 provenance, checked by archive checksum rather than declared as
        # an input, so the known trim-provenance mtime cascade cannot reopen the screens.
        provenance=f"{WORK}/provenance/phyloflash/{{sample}}.json",
    log:
        f"{WORK}/logs/eukaryote_gate/extract/{{sample}}.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --sample {wildcards.sample} --provenance {params.provenance:q} \
          --r1 {output.r1:q} --r2 {output.r2:q} --assembled {output.assembled:q} \
          --classification {output.classification:q} > {log:q} 2>&1
        """


rule build_eukaryote_ssu_reference:
    input:
        script="workflow/scripts/build_eukaryote_ssu_reference.py",
        silva=EUK_SILVA,
        classifications=expand(f"{EUK}/ssu/{{sample}}.ssu_classification.csv", sample=SAMPLE_IDS),
        assembled=expand(f"{EUK}/ssu/{{sample}}.assembled_ssu.fasta", sample=SAMPLE_IDS),
    output:
        fasta=f"{EUK}/reference/ssu_reference.fasta",
        taxonomy=f"{EUK}/reference/ssu_reference.taxonomy.tsv",
    log:
        f"{WORK}/logs/eukaryote_gate/build_eukaryote_ssu_reference.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --silva {input.silva:q} --classifications {input.classifications:q} \
          --fasta {output.fasta:q} --taxonomy {output.taxonomy:q} > {log:q} 2>&1
        """


rule index_eukaryote_ssu_reference_minimap2:
    input:
        f"{EUK}/reference/ssu_reference.fasta",
    output:
        f"{EUK}/reference/ssu_reference.sr.mmi",
    params:
        tools=EUK_TOOLS,
    log:
        f"{WORK}/logs/eukaryote_gate/index_ssu_reference.log",
    threads: 8
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        minimap2 -x sr -t {threads} -d {output:q}.tmp {input:q} > {log:q} 2>&1
        mv {output:q}.tmp {output:q}
        """


rule index_eukaryote_ssu_reference_asm20:
    # A minimap2 index fixes its preset, so full-length sequences need their own index.
    input:
        f"{EUK}/reference/ssu_reference.fasta",
    output:
        f"{EUK}/reference/ssu_reference.asm20.mmi",
    params:
        tools=EUK_TOOLS,
    log:
        f"{WORK}/logs/eukaryote_gate/index_ssu_reference_asm20.log",
    threads: 8
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        minimap2 -x asm20 -t {threads} -d {output:q}.tmp {input:q} > {log:q} 2>&1
        mv {output:q}.tmp {output:q}
        """


rule classify_ssu_reads_minimap2:
    input:
        r1=f"{EUK}/ssu/{{sample}}.R1.fq.gz",
        r2=f"{EUK}/ssu/{{sample}}.R2.fq.gz",
        index=f"{EUK}/reference/ssu_reference.sr.mmi",
        taxonomy=f"{EUK}/reference/ssu_reference.taxonomy.tsv",
        script="workflow/scripts/classify_ssu_reads_lca.py",
    output:
        f"{EUK}/read_lineages/{{sample}}.tsv",
    params:
        tools=EUK_TOOLS,
    log:
        f"{WORK}/logs/eukaryote_gate/classify/{{sample}}.log",
    threads: 8
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        # Keep near-tied secondary hits so each pair is assigned to the LCA of its best hits.
        minimap2 -ax sr -N 200 -p 0.99 --secondary=yes -t {threads} {input.index:q} {input.r1:q} {input.r2:q} \
            2> {log:q} \
          | python {input.script:q} --taxonomy {input.taxonomy:q} --sample {wildcards.sample} \
              --output {output:q} 2>> {log:q}
        """


rule aggregate_ssu_read_lineages:
    input:
        tables=expand(f"{EUK}/read_lineages/{{sample}}.tsv", sample=SAMPLE_IDS),
    output:
        f"{EUK}/read_lineages.tsv",
    log:
        f"{WORK}/logs/eukaryote_gate/aggregate_read_lineages.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        # Keep the header once, then every library's rows.
        awk 'FNR > 1 || NR == 1' {input.tables:q} > {output:q}.tmp
        mv {output:q}.tmp {output:q}
        echo "rows=$(($(wc -l < {output:q}) - 1))" > {log:q}
        """


rule eukaryote_evidence:
    input:
        f"{EUK}/read_lineages.tsv",


rule classify_assembled_ssu_minimap2:
    # One job per library: phyloFlash SSUs, plus 18S genes from the library's Phase-3 assembly
    # when it was assembled, each LCA-classified against the competitive reference.
    input:
        phyloflash=f"{EUK}/ssu/{{sample}}.assembled_ssu.fasta",
        assembly=lambda wildcards: [f"{WORK}/phase3_cohort/markers/{wildcards.sample}.gff",
                                    f"{WORK}/phase3_cohort/assemblies/{wildcards.sample}.fasta"]
                                   if wildcards.sample in ASSEMBLY_IDS else [],
        index=f"{EUK}/reference/ssu_reference.asm20.mmi",
        taxonomy=f"{EUK}/reference/ssu_reference.taxonomy.tsv",
        extract="workflow/scripts/extract_assembly_18s.py",
        classify="workflow/scripts/classify_sequences_lca.py",
    output:
        phyloflash=f"{EUK}/assembled/{{sample}}.phyloflash.tsv",
        assembly18s=f"{EUK}/assembled/{{sample}}.assembly18s.tsv",
        fasta18s=f"{EUK}/assembled/{{sample}}.assembly18s.fasta",
    params:
        tools=EUK_TOOLS,
        min_length=EUKARYOTE["assembled_classification"]["assembly_18s_min_length"],
        # Empty for the libraries below the assembly floor, which have no Phase-3 assembly.
        gff=lambda wildcards, input: input.assembly[0] if input.assembly else "",
        contigs=lambda wildcards, input: input.assembly[1] if input.assembly else "",
    log:
        f"{WORK}/logs/eukaryote_gate/assembled/{{sample}}.log",
    threads: 4
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        : > {log:q}
        if [ -n "{params.gff}" ]; then
          python {input.extract:q} --gff {params.gff:q} --contigs {params.contigs:q} \
            --min-length {params.min_length} --output {output.fasta18s:q} >> {log:q} 2>&1
        else
          : > {output.fasta18s:q}
        fi
        for source in phyloflash assembly18s; do
          if [ "$source" = phyloflash ]; then query={input.phyloflash:q}; out={output.phyloflash:q};
          else query={output.fasta18s:q}; out={output.assembly18s:q}; fi
          if [ -s "$query" ]; then
            minimap2 -ax asm20 -N 200 -p 0.99 --secondary=yes -t {threads} {input.index:q} "$query" 2>> {log:q} \
              | python {input.classify:q} --taxonomy {input.taxonomy:q} --sample {wildcards.sample} \
                  --output "$out" 2>> {log:q}
          else
            printf 'sample_id\tsequence_id\tlineage\tidentity\taligned_bases\n' > "$out"
          fi
        done
        """


rule collect_eukaryote_assembled:
    input:
        script="workflow/scripts/collect_eukaryote_assembled.py",
        tables=[f"{EUK}/assembled/{s}.{src}.tsv" for s in SAMPLE_IDS for src in ("phyloflash", "assembly18s")],
        fastas=[f"{EUK}/ssu/{s}.assembled_ssu.fasta" if src == "phyloflash" else f"{EUK}/assembled/{s}.assembly18s.fasta"
                for s in SAMPLE_IDS for src in ("phyloflash", "assembly18s")],
    output:
        table=f"{EUK}/assembled_nonhost.tsv",
        fasta=f"{EUK}/assembled_nonhost.fasta",
    log:
        f"{WORK}/logs/eukaryote_gate/collect_assembled.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --classifications {input.tables:q} --fastas {input.fastas:q} \
          --table {output.table:q} --fasta {output.fasta:q} > {log:q} 2>&1
        """


rule verify_eukaryotes_ncbi:
    input:
        script="workflow/scripts/verify_eukaryotes_ncbi.py",
        fasta=f"{EUK}/assembled_nonhost.fasta",
    output:
        f"{EUK}/ncbi_verification.tsv",
    params:
        email="casey.dunn@yale.edu",
    log:
        f"{WORK}/logs/eukaryote_gate/verify_ncbi.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --fasta {input.fasta:q} --email {params.email} --output {output:q} > {log:q} 2>&1
        """


rule grade_eukaryotes:
    input:
        script="workflow/scripts/grade_eukaryotes.py",
        reporting="workflow/scripts/eukaryote_reporting.py",
        units="workflow/scripts/eukaryote_units.py",
        config=str(EUKARYOTE_CONFIG_PATH),
        reads=f"{EUK}/read_lineages.tsv",
        assembled=f"{EUK}/assembled_nonhost.tsv",
        verification=f"{EUK}/ncbi_verification.tsv",
    output:
        f"{EUK}/eukaryote_grades.tsv",
    log:
        f"{WORK}/logs/eukaryote_gate/grade.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --config {input.config:q} --reads {input.reads:q} --assembled {input.assembled:q} \
          --verification {input.verification:q} --output {output:q} > {log:q} 2>&1
        """


rule validate_eukaryote_gate:
    input:
        script="workflow/scripts/validate_eukaryote_gate.py",
        config=str(EUKARYOTE_CONFIG_PATH),
        grades=f"{EUK}/eukaryote_grades.tsv",
        assembled=f"{EUK}/assembled_nonhost.tsv",
        verification=f"{EUK}/ncbi_verification.tsv",
    output:
        f"{WORK}/stages/eukaryote_gate.done",
    log:
        f"{WORK}/logs/eukaryote_gate/validate.log",
    threads: 1
    resources:
        **phase5_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --config {input.config:q} --grades {input.grades:q} \
          --assembled {input.assembled:q} --verification {input.verification:q} --output {output:q} > {log:q} 2>&1
        """


rule eukaryote_gate:
    input:
        f"{WORK}/stages/eukaryote_gate.done",
