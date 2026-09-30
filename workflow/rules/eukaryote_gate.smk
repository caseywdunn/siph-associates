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

Next: presence rules calibrated against negative-control lineages grade these
counts together with the assembled SSUs.
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
