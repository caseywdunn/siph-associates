"""Compare protein competition on the same held-out controls; no read rescanning."""

rule prepare_coi_competitor_proteins:
    input:
        f"{RUN}/panel/target.faa", f"{SOURCE}/annotated.faa",
        f"{RUN}/panel/reference_metadata.tsv", "scripts/crustacean_coi/protein_validation.py",
    output:
        f"{RUN}/protein_competition/reference.faa", f"{RUN}/protein_competition/reference_summary.json",
    resources:
        mem_mb=2000,
    shell:
        "{PY:q} scripts/crustacean_coi/protein_validation.py prepare --root {RUN:q} --source {SOURCE:q}"

rule index_coi_competitor_proteins:
    input:
        f"{RUN}/protein_competition/reference.faa",
    output:
        f"{RUN}/protein_competition/reference.dmnd",
    threads: 4
    resources:
        mem_mb=2000,
    shell:
        "{DIAMOND:q} makedb --in {input:q} --db {output:q} --threads {threads}"

rule align_heldout_coi_proteins:
    input:
        query=f"{RUN}/heldout/candidates.fna", db=f"{RUN}/protein_competition/reference.dmnd",
    output:
        f"{RUN}/protein_competition/hits.tsv",
    log:
        f"{RUN}/protein_competition/diamond.log",
    threads: 4
    resources:
        mem_mb=5000,
    shell:
        "{DIAMOND:q} blastx --query {input.query:q} --db {input.db:q} --out {output:q} "
        "--query-gencode 5 --sensitive --min-orf 20 --evalue 1e-3 --query-cover 60 "
        "--max-target-seqs 200 --block-size 0.2 --index-chunks 1 --threads {threads} "
        "--outfmt 6 qseqid sseqid pident length bitscore evalue 2> {log:q}"

rule report_protein_competition:
    input:
        f"{RUN}/protein_competition/hits.tsv", f"{RUN}/heldout/candidate_pairs.tsv",
        f"{RUN}/heldout/truth.json", f"{RUN}/panel/reference_metadata.tsv",
        "scripts/crustacean_coi/protein_validation.py", "scripts/crustacean_coi/summarize.py",
    output:
        f"{RUN}/protein_competition/summary.json", f"{RUN}/protein_competition/classified_controls.tsv",
    shell:
        "{PY:q} scripts/crustacean_coi/protein_validation.py report --root {RUN:q}"

rule protein_validation:
    input:
        rules.report_protein_competition.output,

rule align_all_near_best_coi_proteins:
    input:
        query=f"{RUN}/heldout/candidates.fna", db=f"{RUN}/protein_competition/reference.dmnd",
    output:
        f"{RUN}/protein_competition/all_near_best_hits.tsv",
    log:
        f"{RUN}/protein_competition/all_near_best_diamond.log",
    threads: 4
    resources:
        mem_mb=5000,
    shell:
        "{DIAMOND:q} blastx --query {input.query:q} --db {input.db:q} --out {output:q} "
        "--query-gencode 5 --sensitive --min-orf 20 --evalue 1e-3 --query-cover 60 "
        "--top 1 --block-size 0.2 --index-chunks 1 --threads {threads} "
        "--outfmt 6 qseqid sseqid pident length bitscore evalue 2> {log:q}"

rule report_all_near_best_protein_competition:
    input:
        f"{RUN}/protein_competition/all_near_best_hits.tsv", f"{RUN}/heldout/candidate_pairs.tsv",
        f"{RUN}/heldout/truth.json", f"{RUN}/panel/reference_metadata.tsv",
        "scripts/crustacean_coi/protein_validation.py", "scripts/crustacean_coi/summarize.py",
    output:
        f"{RUN}/protein_competition/all_near_best_summary.json",
        f"{RUN}/protein_competition/all_near_best_classified_controls.tsv",
    shell:
        "{PY:q} scripts/crustacean_coi/protein_validation.py report --root {RUN:q} --all-near-best"

rule protein_validation_all_near_best:
    input:
        rules.report_all_near_best_protein_competition.output,
