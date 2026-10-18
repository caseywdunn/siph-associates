"""Shared competitive protein search for controls, pilots and cohort candidates.

Both genetic codes compete by bit score in the downstream origin report. Retain
all reported hits within 5% of the best score rather than truncating near ties.
"""

rule verify_candidate_origins_diamond:
    input:
        query=coi_verify_query,
        db=coi_verify_db,
        gate=coi_verify_gate,
    output:
        f"{VERIFY_ROOT}/{{sample}}/code{{code}}.tsv",
    wildcard_constraints:
        code="[15]",
    log:
        f"{VERIFY_ROOT}/{{sample}}/code{{code}}.log",
    threads: 4
    resources:
        mem_mb=4000,
    benchmark:
        f"{VERIFY_ROOT}/{{sample}}/code{{code}}.benchmark.tsv",
    shell:
        "if test -s {input.query:q}; then "
        "{DIAMOND:q} blastx --query {input.query:q} --db {input.db:q} --out {output:q} "
        "--query-gencode {wildcards.code} --very-sensitive --min-orf 20 --evalue 1e-3 "
        "--query-cover 60 --top 5 --max-hsps 1 --block-size 0.2 --index-chunks 1 --threads {threads} "
        "--outfmt 6 qseqid sseqid pident length qcovhsp bitscore evalue 2> {log:q}; "
        "else touch {output:q} {log:q}; fi"
