"""
Build the cohort MAG and vOTU catalogs from the accepted per-library assembly
products, applying the locked rules in config/phase3_catalog.json (reasoning in
docs/phase3_qc_decisions.md). Catalog members keep their source-library IDs.

1. Select MAGs by CheckM2 quality: select_mags.
2. Place MAGs in the GTDB r220 taxonomy: prepare_gtdbtk_reference,
   classify_mags_gtdbtk.
3. Keep prokaryotic MAGs and cluster species by ANI: compare_mags_skani,
   build_mag_catalog.
4. Select viral contigs with positive viral evidence: select_viruses.
5. Separate endogenous candidates from associate viruses:
   align_viruses_host_minimap2, flag_endogenous_viruses.
6. Cluster each viral class into vOTUs: cluster_votus_checkv, build_votu_catalog.
7. Check the catalogs against the locked rules: validate_catalog.

Next: Phase 4 pools these representatives with taxon-verified reference
genomes and freezes the mapping catalogs.
"""

CATALOG = f"{WORK}/phase3_catalog"
CATALOG_SCRATCH = f"{SCRATCH}/phase3_catalog"
CATALOG_RESOURCES = PHASE3_CATALOG["resources"]
CATALOG_TOOLS = {name: f"{prefix}/bin" for name, prefix in PHASE3_CATALOG["tool_prefixes"].items()}
VIRUS_SETS = {"associate": "VOTU", "endogenous_candidate": "EVEC"}
HOST_ROUTES = sorted(config["databases"]["host_references"])

wildcard_constraints:
    virus_set="associate|endogenous_candidate",


rule select_mags:
    input:
        script="workflow/scripts/select_mags.py",
        config=str(PHASE3_CATALOG_CONFIG_PATH),
        membership=str(PHASE3_ASSEMBLY_MEMBERSHIP_PATH),
        summary=f"{WORK}/phase3_cohort/assembly_summary.tsv",
    output:
        genomes=directory(f"{CATALOG}/mags/candidates"),
        table=f"{CATALOG}/mags/candidates.tsv",
    params:
        results=lambda wildcards, input: str(Path(input.summary).parent),
    log:
        f"{WORK}/logs/phase3_catalog/select_mags.log",
    threads: 1
    resources:
        **catalog_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --membership {input.membership:q} --results {params.results:q} \
          --config {input.config:q} --genomes {output.genomes:q} --table {output.table:q} > {log:q} 2>&1
        """


rule prepare_gtdbtk_reference:
    # The installed r220 package keeps its reference genomes under skani/genomes/
    # while genome_paths.tsv expects skani/database/; link a corrected view.
    input:
        database=ancient(PHASE3_CATALOG["databases"]["gtdbtk"]),
    output:
        directory(f"{CATALOG_SCRATCH}/gtdbtk_r220"),
    log:
        f"{WORK}/logs/phase3_catalog/prepare_gtdbtk_reference.log",
    threads: 1
    resources:
        **catalog_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        mkdir -p {output:q}/skani
        for part in markers masks metadata mrca_red msa pplacer radii split taxonomy; do
          ln -s {input.database:q}/$part {output:q}/$part
        done
        ln -s {input.database:q}/skani/genome_paths.tsv {output:q}/skani/genome_paths.tsv
        ln -s {input.database:q}/skani/genomes {output:q}/skani/database
        cat {input.database:q}/metadata/metadata.txt > {log:q}
        """


rule classify_mags_gtdbtk:
    input:
        genomes=f"{CATALOG}/mags/candidates",
        # A regenerated scratch view must not trigger reclassification.
        reference=ancient(f"{CATALOG_SCRATCH}/gtdbtk_r220"),
    output:
        bacteria=f"{CATALOG}/mags/gtdbtk.bac120.summary.tsv",
        archaea=f"{CATALOG}/mags/gtdbtk.ar53.summary.tsv",
    params:
        tools=CATALOG_TOOLS["gtdbtk"],
        scratch=f"{CATALOG_SCRATCH}/gtdbtk_tmp",
        mash_db=f"{CATALOG_SCRATCH}/gtdbtk_mash/gtdb_r220.msh",
        pplacer_cpus=8,
    log:
        f"{WORK}/logs/phase3_catalog/classify_mags_gtdbtk.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_catalog/classify_mags_gtdbtk.tsv",
    threads: CATALOG_RESOURCES["gtdbtk"]["threads"]
    resources:
        **catalog_resources("gtdbtk"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH GTDBTK_DATA_PATH={input.reference:q}
        mkdir -p {params.scratch:q} "$(dirname {params.mash_db:q})"
        tmp=$(mktemp -d {params.scratch:q}/run.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        gtdbtk --version > {log:q} 2>&1
        gtdbtk classify_wf --genome_dir {input.genomes:q} --extension fa --out_dir "$tmp/out" \
          --cpus {threads} --pplacer_cpus {params.pplacer_cpus} --mash_db {params.mash_db:q} \
          --tmpdir "$tmp" >> {log:q} 2>&1
        header=$(head -n 1 "$tmp/out/gtdbtk.bac120.summary.tsv" 2>/dev/null || echo user_genome)
        [ -f "$tmp/out/gtdbtk.ar53.summary.tsv" ] || echo "$header" > "$tmp/out/gtdbtk.ar53.summary.tsv"
        [ -f "$tmp/out/gtdbtk.bac120.summary.tsv" ] || echo "$header" > "$tmp/out/gtdbtk.bac120.summary.tsv"
        mv "$tmp/out/gtdbtk.ar53.summary.tsv" {output.archaea:q}
        mv "$tmp/out/gtdbtk.bac120.summary.tsv" {output.bacteria:q}
        """


rule compare_mags_skani:
    input:
        genomes=f"{CATALOG}/mags/candidates",
    output:
        f"{CATALOG}/mags/candidates.skani.tsv",
    params:
        tools=CATALOG_TOOLS["skani"],
    log:
        f"{WORK}/logs/phase3_catalog/compare_mags_skani.log",
    threads: CATALOG_RESOURCES["skani"]["threads"]
    resources:
        **catalog_resources("skani"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        skani --version > {log:q} 2>&1
        skani triangle -t {threads} -E {input.genomes:q}/*.fa > {output:q}.tmp 2>> {log:q}
        mv {output:q}.tmp {output:q}
        """


rule build_mag_catalog:
    input:
        script="workflow/scripts/build_mag_catalog.py",
        config=str(PHASE3_CATALOG_CONFIG_PATH),
        candidates=f"{CATALOG}/mags/candidates.tsv",
        genomes=f"{CATALOG}/mags/candidates",
        bacteria=f"{CATALOG}/mags/gtdbtk.bac120.summary.tsv",
        archaea=f"{CATALOG}/mags/gtdbtk.ar53.summary.tsv",
        ani=f"{CATALOG}/mags/candidates.skani.tsv",
    output:
        catalog=f"{CATALOG}/mags/mag_catalog.tsv",
        representatives=directory(f"{CATALOG}/mags/species_representatives"),
    log:
        f"{WORK}/logs/phase3_catalog/build_mag_catalog.log",
    threads: 1
    resources:
        **catalog_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --candidates {input.candidates:q} --genomes {input.genomes:q} \
          --gtdbtk {input.bacteria:q} {input.archaea:q} --ani {input.ani:q} --config {input.config:q} \
          --catalog {output.catalog:q} --representatives {output.representatives:q} > {log:q} 2>&1
        """


rule select_viruses:
    input:
        script="workflow/scripts/select_viruses.py",
        config=str(PHASE3_CATALOG_CONFIG_PATH),
        membership=str(PHASE3_ASSEMBLY_MEMBERSHIP_PATH),
        summary=f"{WORK}/phase3_cohort/assembly_summary.tsv",
    output:
        fasta=f"{CATALOG}/viruses/included.fna",
        table=f"{CATALOG}/viruses/included.tsv",
    params:
        results=lambda wildcards, input: str(Path(input.summary).parent),
    log:
        f"{WORK}/logs/phase3_catalog/select_viruses.log",
    threads: 1
    resources:
        **catalog_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --membership {input.membership:q} --results {params.results:q} \
          --config {input.config:q} --fasta {output.fasta:q} --table {output.table:q} > {log:q} 2>&1
        """


rule align_viruses_host_minimap2:
    input:
        viruses=f"{CATALOG}/viruses/included.fna",
        reference=lambda wildcards: config["databases"]["host_references"][wildcards.route],
    output:
        f"{CATALOG}/viruses/included.vs_{{route}}.paf",
    params:
        tools=CATALOG_TOOLS["minimap2"],
    log:
        f"{WORK}/logs/phase3_catalog/align_viruses_host_minimap2.{{route}}.log",
    threads: CATALOG_RESOURCES["host_alignment"]["threads"]
    resources:
        **catalog_resources("host_alignment"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        minimap2 -x asm20 -t {threads} --secondary=no -c {input.reference:q} {input.viruses:q} \
          > {output:q}.tmp 2> {log:q}
        mv {output:q}.tmp {output:q}
        """


rule flag_endogenous_viruses:
    input:
        script="workflow/scripts/flag_endogenous_viruses.py",
        config=str(PHASE3_CATALOG_CONFIG_PATH),
        included=f"{CATALOG}/viruses/included.tsv",
        fasta=f"{CATALOG}/viruses/included.fna",
        paf=expand(f"{CATALOG}/viruses/included.vs_{{route}}.paf", route=HOST_ROUTES),
        samples=config["samples"],
        manifest=config["manifest"],
    output:
        table=f"{CATALOG}/viruses/virus_classification.tsv",
        associate=f"{CATALOG}/viruses/associate.fna",
        endogenous=f"{CATALOG}/viruses/endogenous_candidate.fna",
    params:
        paf=lambda wildcards, input: [f"{route}={path}" for route, path in zip(HOST_ROUTES, input.paf)],
    log:
        f"{WORK}/logs/phase3_catalog/flag_endogenous_viruses.log",
    threads: 1
    resources:
        **catalog_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --included {input.included:q} --fasta {input.fasta:q} --paf {params.paf:q} \
          --samples {input.samples:q} --manifest {input.manifest:q} --config {input.config:q} \
          --table {output.table:q} --associate {output.associate:q} --endogenous {output.endogenous:q} \
          > {log:q} 2>&1
        """


rule cluster_votus_checkv:
    input:
        f"{CATALOG}/viruses/{{virus_set}}.fna",
    output:
        ani=f"{CATALOG}/votus/{{virus_set}}.ani.tsv",
        clusters=f"{CATALOG}/votus/{{virus_set}}.clusters.tsv",
    params:
        tools=CATALOG_TOOLS["votu"],
        scratch=f"{CATALOG_SCRATCH}/votu_tmp",
        min_ani=PHASE3_CATALOG["viruses"]["votu_clustering"]["ani_min"],
        min_tcov=PHASE3_CATALOG["viruses"]["votu_clustering"]["shorter_sequence_coverage_min"],
    log:
        f"{WORK}/logs/phase3_catalog/cluster_votus_checkv.{{virus_set}}.log",
    benchmark:
        f"{WORK}/benchmarks/phase3_catalog/cluster_votus_checkv.{{virus_set}}.tsv",
    threads: CATALOG_RESOURCES["votu"]["threads"]
    resources:
        **catalog_resources("votu"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        mkdir -p {params.scratch:q}
        tmp=$(mktemp -d {params.scratch:q}/{wildcards.virus_set}.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        blastn -version > {log:q} 2>&1
        # MIUViG procedure (CheckV): all-vs-all megablast, pairwise ANI, centroid clustering.
        makeblastdb -in {input:q} -dbtype nucl -out "$tmp/db" >> {log:q} 2>&1
        blastn -query {input:q} -db "$tmp/db" -outfmt '6 std qlen slen' -max_target_seqs 10000 \
          -num_threads {threads} -out "$tmp/blast.tsv" 2>> {log:q}
        anicalc -i "$tmp/blast.tsv" -o "$tmp/ani.tsv" >> {log:q} 2>&1
        aniclust --fna {input:q} --ani "$tmp/ani.tsv" --out "$tmp/clusters.tsv" \
          --min_ani {params.min_ani} --min_tcov {params.min_tcov} --min_qcov 0 >> {log:q} 2>&1
        mv "$tmp/ani.tsv" {output.ani:q}
        mv "$tmp/clusters.tsv" {output.clusters:q}
        """


rule build_votu_catalog:
    input:
        script="workflow/scripts/build_votu_catalog.py",
        clusters=f"{CATALOG}/votus/{{virus_set}}.clusters.tsv",
        fasta=f"{CATALOG}/viruses/{{virus_set}}.fna",
        classification=f"{CATALOG}/viruses/virus_classification.tsv",
    output:
        catalog=f"{CATALOG}/votus/{{virus_set}}.votus.tsv",
        representatives=f"{CATALOG}/votus/{{virus_set}}.representatives.fna",
    params:
        prefix=lambda wildcards: VIRUS_SETS[wildcards.virus_set],
    log:
        f"{WORK}/logs/phase3_catalog/build_votu_catalog.{{virus_set}}.log",
    threads: 1
    resources:
        **catalog_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --clusters {input.clusters:q} --fasta {input.fasta:q} \
          --classification {input.classification:q} --prefix {params.prefix} \
          --catalog {output.catalog:q} --representatives {output.representatives:q} > {log:q} 2>&1
        """


rule validate_catalog:
    input:
        script="workflow/scripts/validate_catalog.py",
        config=str(PHASE3_CATALOG_CONFIG_PATH),
        candidates=f"{CATALOG}/mags/candidates.tsv",
        mag_catalog=f"{CATALOG}/mags/mag_catalog.tsv",
        mag_representatives=f"{CATALOG}/mags/species_representatives",
        included=f"{CATALOG}/viruses/included.tsv",
        classification=f"{CATALOG}/viruses/virus_classification.tsv",
        votus=[f"{CATALOG}/votus/{virus_set}.votus.tsv" for virus_set in VIRUS_SETS],
    output:
        f"{WORK}/stages/phase3_catalog.done",
    log:
        f"{WORK}/logs/phase3_catalog/validate_catalog.log",
    threads: 1
    resources:
        **catalog_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --config {input.config:q} --candidates {input.candidates:q} \
          --mag-catalog {input.mag_catalog:q} --mag-representatives {input.mag_representatives:q} \
          --included {input.included:q} --classification {input.classification:q} --votus {input.votus:q} \
          --output {output:q} > {log:q} 2>&1
        """


rule phase3_catalog:
    input:
        f"{WORK}/stages/phase3_catalog.done",
