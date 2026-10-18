"""Extend the tested fish filters to the remaining 156 host-depleted libraries.

Use the pilot implementations unchanged, retaining both mates of their union.
The five completed fish pilots remain external inputs, without repeating scans.
Every new union receives the same competitive whole-mitochondrial BLAST search.
"""

import csv

configfile: "config/fish_12s_host_depleted_cohort.json"

ROOT, PY = config["run"], config["python"]
with open(config["manifest"]) as handle:
    MANIFEST = list(csv.DictReader(handle, delimiter="\t"))
SAMPLES = {r["analysis_sample_id"]: dict(r, index=i) for i, r in enumerate(MANIFEST)}
REUSE = set()
if config.get("reuse_manifest"):
    with open(config["reuse_manifest"]) as handle:
        REUSE = {r["analysis_sample_id"] for r in csv.DictReader(handle, delimiter="\t")}
NEW = sorted(SAMPLES.keys() - REUSE)

wildcard_constraints:
    sample="[^/]+",


rule screen_cohort_fish_seeds:
    input:
        r1=lambda w: SAMPLES[w.sample]["r1_path"],
        r2=lambda w: SAMPLES[w.sample]["r2_path"],
        reference=config["reference"], manifest=config["manifest"],
        code=["scripts/fish_12s/screen_positive.py", "scripts/fish_12s/pilot.py"],
    output:
        directory(f"{ROOT}/screen/{{sample}}/host_depleted"),
    params:
        index=lambda w: SAMPLES[w.sample]["index"],
    threads: 6
    resources:
        mem_mb=8000,
    shell:
        "{PY:q} scripts/fish_12s/screen_positive.py --manifest {config[manifest]:q} "
        "--index {params.index} --threads 4 --reference {input.reference:q} --output-root {ROOT:q}/screen"


rule screen_cohort_fish_aligners:
    input:
        r1=lambda w: SAMPLES[w.sample]["r1_path"],
        r2=lambda w: SAMPLES[w.sample]["r2_path"],
        reference=config["reference"], manifest=config["manifest"],
        index=config["reference_root"] + "/fish_12s.1.bt2",
        code=["scripts/fish_12s/compare_aligners.py", "scripts/fish_12s/pilot.py"],
    output:
        directory(f"{ROOT}/aligners/{{sample}}/host_depleted"),
    params:
        index=lambda w: SAMPLES[w.sample]["index"],
    threads: 6
    resources:
        mem_mb=12000,
    shell:
        "{PY:q} scripts/fish_12s/compare_aligners.py --manifest {config[manifest]:q} "
        "--index {params.index} --threads 4 --reference {input.reference:q} "
        "--bowtie-index {config[reference_root]:q}/fish_12s --sortmerna {config[sortmerna]:q} "
        "--output-root {ROOT:q}/aligners"


rule collect_cohort_fish_pairs:
    input:
        rules.screen_cohort_fish_seeds.output,
        rules.screen_cohort_fish_aligners.output,
        config["manifest"], "scripts/fish_12s/collect_cohort_candidates.py", "scripts/fish_12s/pilot.py",
    output:
        f"{ROOT}/union/{{sample}}/host_depleted/candidates.fasta",
        f"{ROOT}/union/{{sample}}/host_depleted/candidate_pairs.tsv",
        f"{ROOT}/union/{{sample}}/host_depleted/run.json",
    resources:
        mem_mb=2000,
    shell:
        "{PY:q} scripts/fish_12s/collect_cohort_candidates.py --manifest {config[manifest]:q} "
        "--sample {wildcards.sample:q} --screen {ROOT:q}/screen/{wildcards.sample:q}/host_depleted "
        "--aligners {ROOT:q}/aligners/{wildcards.sample:q}/host_depleted "
        "--output {ROOT:q}/union/{wildcards.sample:q}/host_depleted"


rule align_cohort_fish_candidates_blastn:
    input:
        query=f"{ROOT}/union/{{sample}}/host_depleted/candidates.fasta",
        db=config["reference_root"] + "/competitive.nsq",
    output:
        f"{ROOT}/union/{{sample}}/host_depleted/alignment/hits.tsv",
    log:
        f"{ROOT}/union/{{sample}}/host_depleted/alignment/blast.log",
    threads: 4
    resources:
        mem_mb=4000,
    shell:
        "if test -s {input.query:q}; then blastn -task blastn -word_size 11 -dust yes -soft_masking true "
        "-query {input.query:q} -db {config[reference_root]:q}/competitive -evalue 1e-8 "
        "-max_target_seqs 200 -max_hsps 1 -num_threads {threads} "
        "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle' "
        "-out {output:q} 2> {log:q}; else touch {output:q} {log:q}; fi"


rule validate_fish_cohort_controls:
    input:
        f"{ROOT}/union/synthetic/host_depleted/alignment/hits.tsv",
        f"{ROOT}/union/synthetic/host_depleted/candidates.fasta",
        f"{ROOT}/union/synthetic/host_depleted/candidate_pairs.tsv",
        f"{ROOT}/union/synthetic/host_depleted/run.json",
        truth="data/results/fish_12s/positive_aligner_controls/2026-10-17/synthetic/known_answer/checks.json",
        code="scripts/fish_12s/check_cohort_controls.py",
    output:
        f"{ROOT}/CONTROLS_PASSED",
    shell:
        "{PY:q} {input.code:q} --root {ROOT:q} --truth {input.truth:q}"


rule host_depleted:
    input:
        expand(f"{ROOT}/union/{{sample}}/host_depleted/alignment/hits.tsv", sample=NEW),
        expand(config["reuse_screen"] + "/{sample}/host_depleted/alignment/hits.tsv", sample=REUSE),
        expand(config["reuse_aligners"] + "/{sample}/host_depleted/run.json", sample=REUSE),
        ["data/results/fish_12s/host_depleted_aligner_only/2026-10-17/Church2025__YPM-IZ-110907/host_depleted/alignment/hits.tsv"] if REUSE else [],
        [f"{ROOT}/CONTROLS_PASSED"] if config.get("synthetic_validation") else [],
    default_target: True
