"""Audit MAGSP0026 and its related bins against all domains of life.

1. Validate and preserve bins; tile their DNA: prepare_audit_inputs.
2. Obtain Swiss-Prot and taxonomy: retrieve_audit_references, prepare_audit_references.
3. Recover predicted proteins and original marker families: identify_mag_markers_gtdbtk,
   collect_audit_genes.
4. Check known reference proteins: search_audit_controls_diamond, validate_audit_controls.
5. Search proteins and translated DNA: search_mag_genes_swissprot_diamond,
   search_mag_tiles_swissprot_diamond, search_mag_genes_nr_diamond.
6. Report gene, contig and marker support: report_audit_assignments.

Outputs are evidence for manual reassessment; no catalog or manuscript is modified.
"""

configfile: "config/mag_domain_audit.json"

ROOT = config["run"]
PYTHON = config["python"]
DIAMOND = config["diamond"]
CONFIG = config.get("settings_path", "config/mag_domain_audit.json")
SCRIPT = "scripts/mag_domain_audit"
REF = config.get("reference_root", f"{ROOT}/reference")
OUTPUT_FIELDS = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen stitle"
REPORTS = ["genes_swissprot", "tiles_swissprot", "genes_nr"]


rule prepare_audit_inputs:
    input:
        script=f"{SCRIPT}/prepare.py",
        settings=CONFIG,
        bins=[row["fasta"] for row in config["bins"]],
    output:
        bins=directory(f"{ROOT}/inputs"),
        tiles=f"{ROOT}/tiles.fna",
        metadata=f"{ROOT}/tiles.tsv",
        contigs=f"{ROOT}/contigs.tsv",
        provenance=f"{ROOT}/inputs.json",
    resources: mem_mb=2000
    shell:
        "{PYTHON:q} {input.script:q} inputs --config {input.settings:q}"


rule retrieve_audit_references:
    input: settings=CONFIG
    output:
        fasta=f"{REF}/uniprot_sprot.fasta.gz",
        release=f"{REF}/uniprot_release.txt",
        taxonomy=directory(f"{REF}/taxonomy"),
        hashes=f"{REF}/download_sha256.txt",
    params:
        fasta=config["swissprot_url"],
        release=config["swissprot_release_url"],
        taxonomy=config["taxonomy_url"],
    resources: mem_mb=2000
    shell:
        r"""
        mkdir -p {REF:q}/taxonomy
        curl --fail --location --retry 3 -D {REF:q}/swissprot.headers.txt {params.fasta:q} -o {output.fasta:q}.part
        gzip -t {output.fasta:q}.part
        mv {output.fasta:q}.part {output.fasta:q}
        curl --fail --location --retry 3 {params.release:q} -o {output.release:q}
        curl --fail --location --retry 3 -D {REF:q}/taxonomy.headers.txt {params.taxonomy:q} -o {REF:q}/taxdump.tar.gz.part
        gzip -t {REF:q}/taxdump.tar.gz.part
        mv {REF:q}/taxdump.tar.gz.part {REF:q}/taxdump.tar.gz
        tar -xzf {REF:q}/taxdump.tar.gz -C {output.taxonomy:q} nodes.dmp names.dmp merged.dmp
        sha256sum {output.fasta:q} {output.release:q} {REF:q}/taxdump.tar.gz {output.taxonomy:q}/*.dmp > {output.hashes:q}
        """


rule prepare_audit_references:
    input:
        script=f"{SCRIPT}/prepare.py",
        settings=CONFIG,
        fasta=f"{REF}/uniprot_sprot.fasta.gz",
        taxonomy=f"{REF}/taxonomy",
    output:
        metadata=f"{REF}/swissprot_metadata.tsv",
        controls=f"{REF}/controls.faa",
        control_table=f"{REF}/controls.tsv",
        database=f"{REF}/swissprot.dmnd",
    threads: 4
    resources: mem_mb=8000
    shell:
        r"""
        {PYTHON:q} {input.script:q} reference --config {input.settings:q}
        {DIAMOND:q} makedb --in {input.fasta:q} --db {output.database:q} --threads {threads}
        """


rule identify_mag_markers_gtdbtk:
    input:
        bins=f"{ROOT}/inputs",
        database=ancient(config["gtdbtk_data"]),
    output: directory(f"{ROOT}/gtdb_identify")
    params: tools=config["gtdbtk_prefix"] + "/bin"
    log: f"{ROOT}/logs/gtdb_identify.log"
    threads: 8
    resources: mem_mb=16000
    shell:
        r"""
        export PATH={params.tools:q}:$PATH GTDBTK_DATA_PATH={input.database:q}
        gtdbtk identify --genome_dir {input.bins:q} --extension fa --out_dir {output:q} --cpus {threads} > {log:q} 2>&1
        """


rule collect_audit_genes:
    input:
        script=f"{SCRIPT}/prepare.py",
        settings=CONFIG,
        identify=f"{ROOT}/gtdb_identify",
    output:
        proteins=f"{ROOT}/genes.faa",
        metadata=f"{ROOT}/genes.tsv",
        markers=f"{ROOT}/archaeal_marker_hits.tsv",
    resources: mem_mb=2000
    shell:
        "{PYTHON:q} {input.script:q} genes --config {input.settings:q}"


rule search_audit_controls_diamond:
    input:
        query=f"{REF}/controls.faa",
        database=f"{REF}/swissprot.dmnd",
    output: f"{ROOT}/controls/hits.tsv"
    threads: 4
    resources: mem_mb=4000
    shell:
        "{DIAMOND:q} blastp --query {input.query:q} --db {input.database:q} --out {output:q} --threads {threads} --very-sensitive --evalue 1e-5 --top 10 --max-hsps 1 --outfmt 6 {OUTPUT_FIELDS}"


rule validate_audit_controls:
    input:
        script=f"{SCRIPT}/classify.py",
        helper=f"{SCRIPT}/prepare.py",
        settings=CONFIG,
        queries=f"{REF}/controls.tsv",
        hits=f"{ROOT}/controls/hits.tsv",
        metadata=f"{REF}/swissprot_metadata.tsv",
    output:
        table=f"{ROOT}/controls/assignments.tsv",
        passed=f"{ROOT}/controls/CONTROLS_PASSED",
    resources: mem_mb=4000
    shell:
        "{PYTHON:q} {input.script:q} --config {input.settings:q} --queries {input.queries:q} --hits {input.hits:q} --source swissprot --output {output.table:q} --controls"


rule search_mag_genes_swissprot_diamond:
    input:
        query=f"{ROOT}/genes.faa",
        database=f"{REF}/swissprot.dmnd",
        controls=f"{ROOT}/controls/CONTROLS_PASSED",
    output: f"{ROOT}/hits/genes_swissprot.tsv"
    log: f"{ROOT}/logs/genes_swissprot.log"
    priority: 10
    threads: 8
    resources: mem_mb=12000
    shell:
        "{DIAMOND:q} blastp --query {input.query:q} --db {input.database:q} --out {output:q} --threads {threads} --very-sensitive --evalue 1e-5 --top 10 --max-hsps 1 --outfmt 6 {OUTPUT_FIELDS} > {log:q} 2>&1"


rule search_mag_tiles_swissprot_diamond:
    input:
        query=f"{ROOT}/tiles.fna",
        database=f"{REF}/swissprot.dmnd",
        controls=f"{ROOT}/controls/CONTROLS_PASSED",
    output: f"{ROOT}/hits/tiles_swissprot.tsv"
    log: f"{ROOT}/logs/tiles_swissprot.log"
    priority: 10
    threads: 8
    resources: mem_mb=12000
    shell:
        "{DIAMOND:q} blastx --query {input.query:q} --db {input.database:q} --out {output:q} --threads {threads} --very-sensitive --evalue 1e-5 --top 10 --max-hsps 1 --outfmt 6 {OUTPUT_FIELDS} > {log:q} 2>&1"


rule record_nr_reference:
    input: config["nr_database"]
    output: f"{REF}/nr_sha256.txt"
    resources: mem_mb=1000
    shell:
        "sha256sum {input:q} > {output:q}"


rule search_mag_genes_nr_diamond:
    input:
        query=f"{ROOT}/genes.faa",
        database=config["nr_database"],
        reference_hash=f"{REF}/nr_sha256.txt",
        controls=f"{ROOT}/controls/CONTROLS_PASSED",
    output: f"{ROOT}/hits/genes_nr.tsv"
    log: f"{ROOT}/logs/genes_nr.log"
    params:
        sensitivity=config.get("nr_sensitivity", "--very-sensitive"),
        block_size=config.get("nr_block_size", 2),
        index_chunks=config.get("nr_index_chunks", 4),
    threads: 16
    resources: mem_mb=48000
    shell:
        r"""
        {DIAMOND:q} blastp --query {input.query:q} --db {input.database:q} --out {output:q} \
          --threads {threads} {params.sensitivity} --evalue 1e-5 --top 10 --max-hsps 1 \
          --block-size {params.block_size} --index-chunks {params.index_chunks} --tmpdir "${{TMPDIR:-/tmp}}" \
          --outfmt 6 {OUTPUT_FIELDS} > {log:q} 2>&1
        """


rule report_audit_assignments:
    input:
        script=f"{SCRIPT}/classify.py",
        helper=f"{SCRIPT}/prepare.py",
        settings=CONFIG,
        queries=lambda w: f"{ROOT}/{'tiles' if w.review == 'tiles_swissprot' else 'genes'}.tsv",
        hits=f"{ROOT}/hits/{{review}}.tsv",
        markers=f"{ROOT}/archaeal_marker_hits.tsv",
        metadata=f"{REF}/swissprot_metadata.tsv",
        taxonomy=f"{REF}/taxonomy",
    output:
        table=f"{ROOT}/review/{{review}}.tsv",
        counts=f"{ROOT}/review/{{review}}.counts.tsv",
    params: source=lambda w: "nr" if w.review == "genes_nr" else "swissprot"
    wildcard_constraints: review="genes_swissprot|tiles_swissprot|genes_nr"
    priority: 20
    resources: mem_mb=12000
    shell:
        "{PYTHON:q} {input.script:q} --config {input.settings:q} --queries {input.queries:q} --hits {input.hits:q} --source {params.source:q} --output {output.table:q}"


rule all:
    input:
        expand(f"{ROOT}/review/{{review}}.counts.tsv", review=REPORTS),
        f"{ROOT}/controls/CONTROLS_PASSED",
    default_target: True
