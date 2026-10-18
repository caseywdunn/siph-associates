"""Shared COI recruitment and competitive assignment for pilot and validation."""

rule mask_coi_seeds:
    input:
        f"{PANEL}/target.fna",
    output:
        f"{PANEL}/target.masked.fna",
    resources:
        mem_mb=2000,
    shell:
        "{BLAST:q}/dustmasker -in {input:q} -outfmt fasta | "
        "{SEQKIT:q} replace --threads 1 -s -p '[a-z]' -r N > {output:q}"

rule index_target_coi_proteins:
    input:
        f"{PANEL}/target.faa",
    output:
        f"{PANEL}/target.dmnd",
    threads: 4
    resources:
        mem_mb=2000,
    shell:
        "{DIAMOND:q} makedb --in {input:q} --db {output:q} --threads {threads}"

rule index_competitive_coi:
    input:
        f"{PANEL}/competitive.fna",
    output:
        f"{PANEL}/competitive.nsq",
        f"{PANEL}/competitive.nhr",
        f"{PANEL}/competitive.nin",
    shell:
        "{BLAST:q}/makeblastdb -in {input:q} -dbtype nucl -parse_seqids "
        "-blastdb_version 4 -out {PANEL:q}/competitive"


rule nucleotide_coi_recruitment:
    input:
        r1=lambda w: reads(w, 1),
        r2=lambda w: reads(w, 2),
        reference=f"{PANEL}/target.masked.fna",
        controls=control_gate,
        manifest=config["manifest"],
    output:
        r1=f"{RUN}/{{sample}}/bbduk_R1.fastq",
        r2=f"{RUN}/{{sample}}/bbduk_R2.fastq",
        stats=f"{RUN}/{{sample}}/bbduk_stats.tsv",
    log:
        f"{RUN}/{{sample}}/bbduk.log",
    threads: 4
    resources:
        mem_mb=9000,
    benchmark:
        f"{RUN}/{{sample}}/benchmarks/bbduk.tsv"
    shell:
        "{BB:q}/bbduk.sh -Xmx8g in={input.r1:q} in2={input.r2:q} "
        "ref={input.reference:q} outm={output.r1:q} outm2={output.r2:q} "
        "stats={output.stats:q} k=25 hdist=0 maskmiddle=f rcomp=t "
        "rieb=t qtrim=f minlen=0 overwrite=t threads={threads} 2> {log:q}"


rule translated_coi_recruitment:
    input:
        reads=lambda w: reads(w, w.mate),
        db=f"{PANEL}/target.dmnd",
        controls=control_gate,
        manifest=config["manifest"],
    output:
        f"{RUN}/{{sample}}/diamond_R{{mate}}.tsv",
    wildcard_constraints:
        mate="[12]",
    log:
        f"{RUN}/{{sample}}/diamond_R{{mate}}.log",
    threads: 4
    resources:
        mem_mb=5000,
    benchmark:
        f"{RUN}/{{sample}}/benchmarks/diamond_R{{mate}}.tsv"
    shell:
        "{DIAMOND:q} blastx --query {input.reads:q} --db {input.db:q} --out {output:q} "
        "--query-gencode 5 --sensitive --min-orf 20 --evalue 1e-3 --query-cover 60 "
        "--max-target-seqs 1 --block-size 0.2 --index-chunks 1 --threads {threads} "
        "--outfmt 6 qseqid sseqid pident length bitscore evalue 2> {log:q}"


rule collect_coi_candidates:
    input:
        r1=lambda w: reads(w, 1),
        r2=lambda w: reads(w, 2),
        bbduk=f"{RUN}/{{sample}}/bbduk_R1.fastq",
        diamond1=f"{RUN}/{{sample}}/diamond_R1.tsv",
        diamond2=f"{RUN}/{{sample}}/diamond_R2.tsv",
        code="scripts/crustacean_coi/collect_candidates.py",
    output:
        fasta=f"{RUN}/{{sample}}/candidates.fna",
        pairs=f"{RUN}/{{sample}}/candidate_pairs.tsv",
        metrics=f"{RUN}/{{sample}}/recruitment.json",
    params:
        pairs=expected_pairs,
        directory=lambda w: f"{RUN}/{w.sample}",
    resources:
        mem_mb=2000,
    benchmark:
        f"{RUN}/{{sample}}/benchmarks/collect.tsv"
    shell:
        "{PY:q} {input.code:q} --r1 {input.r1:q} --r2 {input.r2:q} "
        "--bbduk {input.bbduk:q} --diamond {input.diamond1:q} {input.diamond2:q} "
        "{params.pairs} --output {params.directory:q}"


rule competitive_coi_alignment:
    input:
        query=f"{RUN}/{{sample}}/candidates.fna",
        db=rules.index_competitive_coi.output,
    output:
        f"{RUN}/{{sample}}/competitive.tsv",
    log:
        f"{RUN}/{{sample}}/competitive.log",
    threads: 4
    resources:
        mem_mb=2000,
    benchmark:
        f"{RUN}/{{sample}}/benchmarks/competitive.tsv"
    shell:
        "if test -s {input.query:q}; then "
        "{BLAST:q}/blastn -query {input.query:q} -db {PANEL:q}/competitive "
        "-task blastn -word_size 11 -dust yes -evalue 1e-8 -max_hsps 1 "
        "-max_target_seqs 200 -num_threads {threads} -out {output:q} "
        "-outfmt '6 qseqid sseqid pident length qstart qend sstart send bitscore evalue' "
        "2> {log:q}; else touch {output:q} {log:q}; fi"


rule summarize_coi:
    input:
        alignment=f"{RUN}/{{sample}}/competitive.tsv",
        pairs=f"{RUN}/{{sample}}/candidate_pairs.tsv",
        metrics=f"{RUN}/{{sample}}/recruitment.json",
        metadata=f"{PANEL}/reference_metadata.tsv",
        code="scripts/crustacean_coi/summarize.py",
    output:
        f"{RUN}/{{sample}}/classified_pairs.tsv",
        f"{RUN}/{{sample}}/summary.json",
    params:
        directory=lambda w: f"{RUN}/{w.sample}",
    shell:
        "{PY:q} {input.code:q} --directory {params.directory:q} --metadata {input.metadata:q}"
