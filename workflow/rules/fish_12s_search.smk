"""Shared fish filters, paired candidate recovery, and competitive alignment."""

rule screen_cohort_fish_seeds:
    input:
        r1=lambda w: fish_read(w, 1),
        r2=lambda w: fish_read(w, 2),
        reference=config["reference"], manifest=fish_manifest,
        code=["scripts/fish_12s/screen_positive.py", "scripts/fish_12s/pilot.py"],
    output:
        directory(f"{ROOT}/screen/{{sample}}/{READ_LABEL}"),
    params:
        index=lambda w: fish_index(w),
    threads: 6
    resources:
        mem_mb=8000,
    shell:
        "{PY:q} scripts/fish_12s/screen_positive.py --manifest {input.manifest:q} "
        "--index {params.index} --threads 4 --reference {input.reference:q} --output-root {ROOT:q}/screen"


rule screen_cohort_fish_aligners:
    input:
        r1=lambda w: fish_read(w, 1),
        r2=lambda w: fish_read(w, 2),
        reference=config["reference"], manifest=fish_manifest,
        index=config["reference_root"] + "/fish_12s.1.bt2",
        code=["scripts/fish_12s/compare_aligners.py", "scripts/fish_12s/pilot.py"],
    output:
        directory(f"{ROOT}/aligners/{{sample}}/{READ_LABEL}"),
    params:
        index=lambda w: fish_index(w),
    threads: 6
    resources:
        mem_mb=12000,
    shell:
        "{PY:q} scripts/fish_12s/compare_aligners.py --manifest {input.manifest:q} "
        "--index {params.index} --threads 4 --reference {input.reference:q} "
        "--bowtie-index {config[reference_root]:q}/fish_12s --sortmerna {config[sortmerna]:q} "
        "--output-root {ROOT:q}/aligners"


rule collect_cohort_fish_pairs:
    input:
        rules.screen_cohort_fish_seeds.output,
        rules.screen_cohort_fish_aligners.output,
        "scripts/fish_12s/collect_cohort_candidates.py", "scripts/fish_12s/pilot.py",
        manifest=fish_manifest,
        r1=lambda w: fish_read(w, 1),
        r2=lambda w: fish_read(w, 2),
    output:
        f"{ROOT}/union/{{sample}}/{READ_LABEL}/candidates.fasta",
        f"{ROOT}/union/{{sample}}/{READ_LABEL}/candidate_pairs.tsv",
        f"{ROOT}/union/{{sample}}/{READ_LABEL}/run.json",
    resources:
        mem_mb=2000,
    shell:
        "{PY:q} scripts/fish_12s/collect_cohort_candidates.py --manifest {input.manifest:q} "
        "--sample {wildcards.sample:q} --screen {ROOT:q}/screen/{wildcards.sample:q}/{READ_LABEL} "
        "--aligners {ROOT:q}/aligners/{wildcards.sample:q}/{READ_LABEL} "
        "--output {ROOT:q}/union/{wildcards.sample:q}/{READ_LABEL}"


rule align_cohort_fish_candidates_blastn:
    input:
        query=f"{ROOT}/union/{{sample}}/{READ_LABEL}/candidates.fasta",
        db=config["reference_root"] + "/competitive.nsq",
    output:
        f"{ROOT}/union/{{sample}}/{READ_LABEL}/alignment/hits.tsv",
    log:
        f"{ROOT}/union/{{sample}}/{READ_LABEL}/alignment/blast.log",
    threads: 4
    resources:
        mem_mb=4000,
    shell:
        "if test -s {input.query:q}; then blastn -task blastn -word_size 11 -dust yes -soft_masking true "
        "-query {input.query:q} -db {config[reference_root]:q}/competitive -evalue 1e-8 "
        "-max_target_seqs 200 -max_hsps 1 -num_threads {threads} "
        "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle' "
        "-out {output:q} 2> {log:q}; else touch {output:q} {log:q}; fi"
