"""
Build and freeze the version-v1 mapping catalogs, applying the locked rules in
config/phase4_catalog.json (reasoning in docs/phase4_catalog_decisions.md).
The bacterial/archaeal catalog pools sylph-nominated reference genomes with the
Phase-3 MAG species and carries 20 decoy genomes as presence-rule controls.

1. Nominate references from cohort sylph hits and draw decoy candidates:
   nominate_references.
2. Obtain genomes: fetch_oceandna_archive, collect_reference_genomes.
3. Screen references: check_references_ncbi, assess_references_checkm2,
   classify_oceandna_gtdbtk.
4. Pool, dereplicate, add decoys, and freeze the catalog:
   compare_catalog_genomes_skani, build_bacterial_catalog,
   index_bacterial_catalog_bwa.
5. Freeze the viral catalogs: build_viral_catalog.
6. Link vOTUs to hosts by CRISPR spacers: call_crispr_spacers_minced,
   match_spacers_votus_blastn.
7. Check the frozen catalogs: validate_phase4_catalog.

Next: Phase 5 maps every library against these frozen catalogs.
"""

P4 = f"{WORK}/phase4_catalog"
P4_SCRATCH = f"{SCRATCH}/phase4_catalog"
P4_TOOLS = {name: f"{prefix}/bin" for name, prefix in PHASE4_CATALOG["tool_prefixes"].items()}
P4_VERSION = PHASE4_CATALOG["catalog_version"]
P4_BACTERIAL = f"{P4}/{P4_VERSION}/bacterial_catalog"
P4_VIRAL_SETS = ("associate", "endogenous_candidate")

wildcard_constraints:
    viral_set="associate|endogenous_candidate",


rule nominate_references:
    input:
        script="workflow/scripts/nominate_references.py",
        config=str(PHASE4_CATALOG_CONFIG_PATH),
        taxonomy=PHASE4_CATALOG["references"]["gtdb_taxonomy"],
    output:
        references=f"{P4}/references/nominated.tsv",
        decoys=f"{P4}/references/decoy_candidates.tsv",
    params:
        # Accepted Phase-2 output, checked against its frozen checksum. Declaring it
        # as an input would pull in the known trim-provenance mtime cascade.
        nominations=f"{WORK}/aggregation/cohort/candidate_nominations.tsv",
        sha256=PHASE4_CATALOG["nomination"]["source_sha256"],
    log:
        f"{WORK}/logs/phase4_catalog/nominate_references.log",
    threads: 1
    resources:
        **phase4_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        echo "{params.sha256}  {params.nominations}" | sha256sum -c - > {log:q} 2>&1
        python {input.script:q} --nominations {params.nominations:q} --taxonomy {input.taxonomy:q} \
          --config {input.config:q} --references {output.references:q} --decoys {output.decoys:q} >> {log:q} 2>&1
        """


rule fetch_oceandna_archive:
    output:
        f"{P4}/oceandna/fasta_species-representatives.tar",
    params:
        url=PHASE4_CATALOG["references"]["oceandna_species_representatives"]["url"],
        md5=PHASE4_CATALOG["references"]["oceandna_species_representatives"]["md5"],
    log:
        f"{WORK}/logs/phase4_catalog/fetch_oceandna_archive.log",
    threads: 1
    resources:
        **phase4_resources("download"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        curl -fsSL --retry 5 -o {output:q}.part {params.url:q} 2> {log:q}
        echo "{params.md5}  {output}.part" | md5sum -c - >> {log:q} 2>&1
        mv {output:q}.part {output:q}
        """


rule collect_reference_genomes:
    input:
        script="workflow/scripts/collect_reference_genomes.py",
        references=f"{P4}/references/nominated.tsv",
        decoys=f"{P4}/references/decoy_candidates.tsv",
        archive=f"{P4}/oceandna/fasta_species-representatives.tar",
    output:
        genomes=directory(f"{P4}/references/genomes"),
        table=f"{P4}/references/staged.tsv",
        oceandna_batch=f"{P4}/references/oceandna_batch.tsv",
    params:
        gtdb_genomes=PHASE4_CATALOG["references"]["gtdb_genomes"],
    log:
        f"{WORK}/logs/phase4_catalog/collect_reference_genomes.log",
    threads: 1
    resources:
        **phase4_resources("download"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --references {input.references:q} --decoys {input.decoys:q} \
          --gtdb-genomes {params.gtdb_genomes:q} --oceandna-archive {input.archive:q} \
          --genomes {output.genomes:q} --table {output.table:q} --oceandna-batch {output.oceandna_batch:q} \
          > {log:q} 2>&1
        """


rule check_references_ncbi:
    input:
        script="workflow/scripts/check_ncbi_status.py",
        staged=f"{P4}/references/staged.tsv",
    output:
        f"{P4}/references/ncbi_status.tsv",
    params:
        email="casey.dunn@yale.edu",
    log:
        f"{WORK}/logs/phase4_catalog/check_references_ncbi.log",
    threads: 1
    resources:
        **phase4_resources("download"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --genomes {input.staged:q} --email {params.email} --output {output:q} > {log:q} 2>&1
        """


rule assess_references_checkm2:
    input:
        genomes=f"{P4}/references/genomes",
        database=PHASE4_CATALOG["databases"]["checkm2"],
    output:
        f"{P4}/references/checkm2_quality.tsv",
    params:
        tools=P4_TOOLS["checkm2"],
        scratch=f"{P4_SCRATCH}/checkm2_tmp",
    log:
        f"{WORK}/logs/phase4_catalog/assess_references_checkm2.log",
    benchmark:
        f"{WORK}/benchmarks/phase4_catalog/assess_references_checkm2.tsv",
    threads: PHASE4_CATALOG["resources"]["checkm2"]["threads"]
    resources:
        **phase4_resources("checkm2"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        mkdir -p {params.scratch:q}
        tmp=$(mktemp -d {params.scratch:q}/run.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        checkm2 --version > {log:q} 2>&1
        checkm2 predict --threads {threads} --input {input.genomes:q} --extension fa \
          --database_path {input.database:q} --output-directory "$tmp/out" --force >> {log:q} 2>&1
        cp "$tmp/out/quality_report.tsv" {output:q}.tmp
        mv {output:q}.tmp {output:q}
        """


rule classify_oceandna_gtdbtk:
    input:
        batch=f"{P4}/references/oceandna_batch.tsv",
        genomes=f"{P4}/references/genomes",
        reference=ancient(f"{CATALOG_SCRATCH}/gtdbtk_r220"),
    output:
        bacteria=f"{P4}/references/oceandna_gtdbtk.bac120.summary.tsv",
        archaea=f"{P4}/references/oceandna_gtdbtk.ar53.summary.tsv",
    params:
        tools=P4_TOOLS["gtdbtk"],
        scratch=f"{P4_SCRATCH}/gtdbtk_tmp",
        mash_db=f"{CATALOG_SCRATCH}/gtdbtk_mash/gtdb_r220.msh",
        pplacer_cpus=8,
    log:
        f"{WORK}/logs/phase4_catalog/classify_oceandna_gtdbtk.log",
    benchmark:
        f"{WORK}/benchmarks/phase4_catalog/classify_oceandna_gtdbtk.tsv",
    threads: PHASE4_CATALOG["resources"]["gtdbtk"]["threads"]
    resources:
        **phase4_resources("gtdbtk"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH GTDBTK_DATA_PATH={input.reference:q}
        mkdir -p {params.scratch:q}
        tmp=$(mktemp -d {params.scratch:q}/run.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        gtdbtk --version > {log:q} 2>&1
        gtdbtk classify_wf --batchfile {input.batch:q} --out_dir "$tmp/out" --cpus {threads} \
          --pplacer_cpus {params.pplacer_cpus} --mash_db {params.mash_db:q} >> {log:q} 2>&1
        header=$(head -n 1 "$tmp/out/gtdbtk.bac120.summary.tsv" 2>/dev/null || echo user_genome)
        [ -e "$tmp/out/gtdbtk.ar53.summary.tsv" ] || echo "$header" > "$tmp/out/gtdbtk.ar53.summary.tsv"
        [ -e "$tmp/out/gtdbtk.bac120.summary.tsv" ] || echo "$header" > "$tmp/out/gtdbtk.bac120.summary.tsv"
        # Top-level summaries are relative symlinks into $tmp; copy their targets.
        cp -L "$tmp/out/gtdbtk.ar53.summary.tsv" {output.archaea:q}.tmp
        cp -L "$tmp/out/gtdbtk.bac120.summary.tsv" {output.bacteria:q}.tmp
        mv {output.archaea:q}.tmp {output.archaea:q}
        mv {output.bacteria:q}.tmp {output.bacteria:q}
        """


rule compare_catalog_genomes_skani:
    input:
        genomes=f"{P4}/references/genomes",
        mags=f"{WORK}/phase3_catalog/mags/species_representatives",
    output:
        f"{P4}/references/catalog_genomes.skani.tsv",
    params:
        tools=P4_TOOLS["skani"],
    log:
        f"{WORK}/logs/phase4_catalog/compare_catalog_genomes_skani.log",
    threads: PHASE4_CATALOG["resources"]["skani"]["threads"]
    resources:
        **phase4_resources("skani"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        skani --version > {log:q} 2>&1
        skani triangle -t {threads} -E {input.genomes:q}/*.fa {input.mags:q}/*.fa > {output:q}.tmp 2>> {log:q}
        mv {output:q}.tmp {output:q}
        """


rule build_bacterial_catalog:
    input:
        script="workflow/scripts/build_bacterial_catalog.py",
        config=str(PHASE4_CATALOG_CONFIG_PATH),
        references=f"{P4}/references/nominated.tsv",
        decoys=f"{P4}/references/decoy_candidates.tsv",
        staged=f"{P4}/references/staged.tsv",
        genomes=f"{P4}/references/genomes",
        ncbi=f"{P4}/references/ncbi_status.tsv",
        checkm2=f"{P4}/references/checkm2_quality.tsv",
        oceandna=[f"{P4}/references/oceandna_gtdbtk.bac120.summary.tsv",
                  f"{P4}/references/oceandna_gtdbtk.ar53.summary.tsv"],
        mag_catalog=f"{WORK}/phase3_catalog/mags/mag_catalog.tsv",
        mags=f"{WORK}/phase3_catalog/mags/species_representatives",
        ani=f"{P4}/references/catalog_genomes.skani.tsv",
    output:
        manifest=f"{P4_BACTERIAL}.manifest.tsv",
        fasta=f"{P4_BACTERIAL}.fna",
    log:
        f"{WORK}/logs/phase4_catalog/build_bacterial_catalog.log",
    threads: 1
    resources:
        **phase4_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --config {input.config:q} --references {input.references:q} \
          --decoys {input.decoys:q} --staged {input.staged:q} --genomes {input.genomes:q} --ncbi {input.ncbi:q} \
          --checkm2 {input.checkm2:q} --oceandna-gtdbtk {input.oceandna:q} --mag-catalog {input.mag_catalog:q} \
          --mag-representatives {input.mags:q} --ani {input.ani:q} --manifest {output.manifest:q} \
          --fasta {output.fasta:q} > {log:q} 2>&1
        """


rule index_bacterial_catalog_bwa:
    input:
        f"{P4_BACTERIAL}.fna",
    output:
        multiext(f"{P4_BACTERIAL}.fna", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"),
    params:
        tools=P4_TOOLS["bwa"],
    log:
        f"{WORK}/logs/phase4_catalog/index_bacterial_catalog_bwa.log",
    benchmark:
        f"{WORK}/benchmarks/phase4_catalog/index_bacterial_catalog_bwa.tsv",
    threads: PHASE4_CATALOG["resources"]["index"]["threads"]
    resources:
        **phase4_resources("index"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        bwa index {input:q} > {log:q} 2>&1
        samtools faidx {input:q} >> {log:q} 2>&1
        """


rule build_viral_catalog:
    input:
        script="workflow/scripts/build_viral_catalog.py",
        representatives=f"{WORK}/phase3_catalog/votus/{{viral_set}}.representatives.fna",
        votus=f"{WORK}/phase3_catalog/votus/{{viral_set}}.votus.tsv",
    output:
        fasta=f"{P4}/{P4_VERSION}/viral_{{viral_set}}.fna",
        manifest=f"{P4}/{P4_VERSION}/viral_{{viral_set}}.manifest.tsv",
    log:
        f"{WORK}/logs/phase4_catalog/build_viral_catalog.{{viral_set}}.log",
    threads: 1
    resources:
        **phase4_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --representatives {input.representatives:q} --votus {input.votus:q} \
          --catalog-class {wildcards.viral_set} --fasta {output.fasta:q} --manifest {output.manifest:q} > {log:q} 2>&1
        """


rule call_crispr_spacers_minced:
    input:
        f"{P4_BACTERIAL}.fna",
    output:
        spacers=f"{P4}/{P4_VERSION}/crispr/spacers.fna",
        arrays=f"{P4}/{P4_VERSION}/crispr/arrays.gff",
    params:
        tools=P4_TOOLS["crispr"],
        scratch=f"{P4_SCRATCH}/minced_tmp",
    log:
        f"{WORK}/logs/phase4_catalog/call_crispr_spacers_minced.log",
    threads: 1
    resources:
        **phase4_resources("crispr"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        mkdir -p {params.scratch:q}
        tmp=$(mktemp -d {params.scratch:q}/run.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        minced --version > {log:q} 2>&1
        minced -spacers -gffFull {input:q} "$tmp/arrays.txt" "$tmp/arrays.gff" >> {log:q} 2>&1
        [ -e "$tmp/arrays_spacers.fa" ] || : > "$tmp/arrays_spacers.fa"
        mv "$tmp/arrays_spacers.fa" {output.spacers:q}
        mv "$tmp/arrays.gff" {output.arrays:q}
        """


rule match_spacers_votus_blastn:
    input:
        script="workflow/scripts/link_spacers_to_votus.py",
        config=str(PHASE4_CATALOG_CONFIG_PATH),
        spacers=f"{P4}/{P4_VERSION}/crispr/spacers.fna",
        votus=[f"{P4}/{P4_VERSION}/viral_{s}.fna" for s in P4_VIRAL_SETS],
        manifest=f"{P4_BACTERIAL}.manifest.tsv",
    output:
        matches=f"{P4}/{P4_VERSION}/crispr/spacer_matches.tsv",
        links=f"{P4}/{P4_VERSION}/crispr/host_links.tsv",
    params:
        tools=P4_TOOLS["crispr"],
        scratch=f"{P4_SCRATCH}/spacer_tmp",
    log:
        f"{WORK}/logs/phase4_catalog/match_spacers_votus_blastn.log",
    threads: PHASE4_CATALOG["resources"]["crispr"]["threads"]
    resources:
        **phase4_resources("crispr"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        mkdir -p {params.scratch:q}
        tmp=$(mktemp -d {params.scratch:q}/run.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        blastn -version > {log:q} 2>&1
        cat {input.votus:q} > "$tmp/votus.fna"
        makeblastdb -in "$tmp/votus.fna" -dbtype nucl -out "$tmp/votus" >> {log:q} 2>&1
        : > "$tmp/blast.tsv"
        if [ -s {input.spacers:q} ]; then
          blastn -task blastn-short -query {input.spacers:q} -db "$tmp/votus" -outfmt '6 std qlen slen' \
            -evalue 1 -max_target_seqs 1000 -num_threads {threads} -out "$tmp/blast.tsv" 2>> {log:q}
        fi
        python {input.script:q} --blast "$tmp/blast.tsv" --manifest {input.manifest:q} --config {input.config:q} \
          --matches {output.matches:q} --links {output.links:q} >> {log:q} 2>&1
        """


rule validate_phase4_catalog:
    input:
        script="workflow/scripts/validate_phase4_catalog.py",
        config=str(PHASE4_CATALOG_CONFIG_PATH),
        references=f"{P4}/references/nominated.tsv",
        manifest=f"{P4_BACTERIAL}.manifest.tsv",
        fasta=f"{P4_BACTERIAL}.fna",
        index=multiext(f"{P4_BACTERIAL}.fna", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"),
        mag_catalog=f"{WORK}/phase3_catalog/mags/mag_catalog.tsv",
        viral_manifests=[f"{P4}/{P4_VERSION}/viral_{s}.manifest.tsv" for s in P4_VIRAL_SETS],
        viral_fastas=[f"{P4}/{P4_VERSION}/viral_{s}.fna" for s in P4_VIRAL_SETS],
        links=f"{P4}/{P4_VERSION}/crispr/host_links.tsv",
    output:
        f"{WORK}/stages/phase4_catalog.done",
    log:
        f"{WORK}/logs/phase4_catalog/validate_phase4_catalog.log",
    threads: 1
    resources:
        **phase4_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --config {input.config:q} --references {input.references:q} \
          --manifest {input.manifest:q} --fasta {input.fasta:q} --mag-catalog {input.mag_catalog:q} \
          --viral-manifests {input.viral_manifests:q} --viral-fastas {input.viral_fastas:q} \
          --links {input.links:q} --output {output:q} > {log:q} 2>&1
        """


rule phase4_catalog:
    input:
        f"{WORK}/stages/phase4_catalog.done",
