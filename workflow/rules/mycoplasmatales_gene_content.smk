"""
Predicted gene content of the siphonophore Metamycoplasmataceae compared with
host-associated Mollicutes and all GTDB r220 species representatives of the
family (config/mycoplasmatales_gene_content.json;
docs/mycoplasmatales_gene_content_plan.md).

1. Collect the genome set and screen its quality: collect_gene_content_genomes,
   assess_gene_content_genomes_checkm2.
2. Obtain KOfam profiles with checksums recorded: download_kofam_profiles.
3. Predict genes with genetic code 4: predict_genes_pyrodigal.
4. Assign KEGG orthologs: annotate_proteins_kofamscan (strict and relaxed tiers are
   applied in the summary).
5. Score KEGG modules and tabulate genomes, lineages, and focal functions:
   summarize_gene_content.

Next: interpretation in the manuscript; absence statements use only
near-complete genomes.
"""

GENE = f"{WORK}/mycoplasmatales_gene_content"
GENE_CONFIG_PATH = ROOT / "config" / "mycoplasmatales_gene_content.json"
with GENE_CONFIG_PATH.open() as handle:
    GENE_CONFIG = json.load(handle)
GENE_TOOLS = f"{GENE_CONFIG['tool_prefixes']['gene_content']}/bin"


def gene_resources(name):
    return {key: GENE_CONFIG["resources"][name][key] for key in ("mem_mb", "runtime", "partition")}


rule collect_gene_content_genomes:
    input:
        script="workflow/scripts/collect_gene_content_genomes.py",
        tree_genomes=f"{WORK}/mycoplasmatales_phylogeny/tree_genomes",
        labels=f"{WORK}/mycoplasmatales_phylogeny/tree_labels.tsv",
        reference=ancient(f"{CATALOG_SCRATCH}/gtdbtk_r220"),
    output:
        genomes=directory(f"{GENE}/genomes"),
        table=f"{GENE}/genomes.tsv",
    params:
        family=GENE_CONFIG["family"],
        taxonomy=GENE_CONFIG["gtdb_taxonomy"],
    log:
        f"{WORK}/logs/mycoplasmatales_gene_content/collect.log",
    threads: 1
    resources:
        **gene_resources("collect"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --tree-genomes {input.tree_genomes:q} --labels {input.labels:q} \
          --gtdb {input.reference:q} --taxonomy {params.taxonomy:q} --family {params.family:q} \
          --genomes {output.genomes:q} --table {output.table:q} > {log:q} 2>&1
        """


rule assess_gene_content_genomes_checkm2:
    input:
        genomes=f"{GENE}/genomes",
        database=PHASE4_CATALOG["databases"]["checkm2"],
    output:
        f"{GENE}/checkm2.tsv",
    params:
        tools=P4_TOOLS["checkm2"],
        scratch=f"{SCRATCH}/mycoplasmatales_gene_content/checkm2_tmp",
    log:
        f"{WORK}/logs/mycoplasmatales_gene_content/checkm2.log",
    threads: GENE_CONFIG["resources"]["checkm2"]["threads"]
    resources:
        **gene_resources("checkm2"),
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


rule download_kofam_profiles:
    output:
        profiles=directory(f"{GENE}/kofam/profiles"),
        ko_list=f"{GENE}/kofam/ko_list",
        checksums=f"{GENE}/kofam/download.tsv",
    params:
        profiles_url=GENE_CONFIG["kofam"]["profiles_url"],
        ko_list_url=GENE_CONFIG["kofam"]["ko_list_url"],
        directory=f"{GENE}/kofam",
    log:
        f"{WORK}/logs/mycoplasmatales_gene_content/download_kofam.log",
    threads: 1
    resources:
        **gene_resources("download"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        cd {params.directory:q}
        curl -fsSL -R -o profiles.tar.gz {params.profiles_url:q} 2> {log:q}
        curl -fsSL -R -o ko_list.gz {params.ko_list_url:q} 2>> {log:q}
        printf 'file\turl\tlast_modified\tsha256\n' > download.tsv.tmp
        for pair in "profiles.tar.gz {params.profiles_url}" "ko_list.gz {params.ko_list_url}"; do
          set -- $pair
          printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$(date -u -r "$1" +%Y-%m-%dT%H:%M:%SZ)" "$(sha256sum "$1" | cut -d' ' -f1)" >> download.tsv.tmp
        done
        tar -xzf profiles.tar.gz
        gunzip -c ko_list.gz > ko_list
        rm profiles.tar.gz
        mv download.tsv.tmp download.tsv
        """


rule predict_genes_pyrodigal:
    input:
        genomes=f"{GENE}/genomes",
    output:
        proteins=directory(f"{GENE}/proteins"),
        genes=directory(f"{GENE}/genes"),
    params:
        tools=GENE_TOOLS,
        code=GENE_CONFIG["genetic_code"],
    log:
        f"{WORK}/logs/mycoplasmatales_gene_content/pyrodigal.log",
    threads: GENE_CONFIG["resources"]["genes"]["threads"]
    resources:
        **gene_resources("genes"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        pyrodigal --version > {log:q} 2>&1
        mkdir -p {output.proteins:q} {output.genes:q}
        for genome in {input.genomes:q}/*.fa; do
          name=$(basename "$genome" .fa)
          pyrodigal -i "$genome" -g {params.code} -j {threads} -f gff -o {output.genes:q}/"$name".gff \
            -a {output.proteins:q}/"$name".faa >> {log:q} 2>&1
        done
        """


rule annotate_proteins_kofamscan:
    input:
        proteins=f"{GENE}/proteins",
        profiles=f"{GENE}/kofam/profiles",
        ko_list=f"{GENE}/kofam/ko_list",
    output:
        f"{GENE}/kofamscan.tsv",
    params:
        tools=GENE_TOOLS,
        scratch=f"{SCRATCH}/mycoplasmatales_gene_content/kofamscan_tmp",
    log:
        f"{WORK}/logs/mycoplasmatales_gene_content/kofamscan.log",
    benchmark:
        f"{WORK}/benchmarks/mycoplasmatales_gene_content/kofamscan.tsv",
    threads: GENE_CONFIG["resources"]["kofamscan"]["threads"]
    resources:
        **gene_resources("kofamscan"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        mkdir -p {params.scratch:q}
        tmp=$(mktemp -d {params.scratch:q}/run.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        # Protein IDs carry their genome so one search covers the whole set.
        for faa in {input.proteins:q}/*.faa; do
          name=$(basename "$faa" .faa)
          sed "s/^>/>$name|/" "$faa"
        done | sed 's/\\*$//' > "$tmp/proteins.faa"
        exec_annotation --version > {log:q} 2>&1
        exec_annotation -o "$tmp/hits.tsv" -p {input.profiles:q} -k {input.ko_list:q} --cpu {threads} \
          --tmp-dir "$tmp/work" -f detail-tsv "$tmp/proteins.faa" >> {log:q} 2>&1
        cp "$tmp/hits.tsv" {output:q}.tmp
        mv {output:q}.tmp {output:q}
        """


rule summarize_gene_content:
    input:
        script="workflow/scripts/summarize_gene_content.py",
        table=f"{GENE}/genomes.tsv",
        checkm2=f"{GENE}/checkm2.tsv",
        genomes=f"{GENE}/genomes",
        genes=f"{GENE}/genes",
        kofam=f"{GENE}/kofamscan.tsv",
        ko_list=f"{GENE}/kofam/ko_list",
    output:
        genomes=f"{GENE}/summary/genome_stats.tsv",
        kos=f"{GENE}/summary/ko_matrix.tsv",
        modules=f"{GENE}/summary/module_completeness.tsv",
        lineages=f"{GENE}/summary/lineage_modules.tsv",
        focal=f"{GENE}/summary/focal_functions.tsv",
        # Written by the script beside the relaxed tables, for the strict-threshold sensitivity view.
        kos_strict=f"{GENE}/summary/ko_matrix_strict.tsv",
        modules_strict=f"{GENE}/summary/module_completeness_strict.tsv",
        lineages_strict=f"{GENE}/summary/lineage_modules_strict.tsv",
    params:
        config=str(GENE_CONFIG_PATH),
        tools=GENE_TOOLS,
        scratch=f"{SCRATCH}/mycoplasmatales_gene_content/modules_tmp",
    log:
        f"{WORK}/logs/mycoplasmatales_gene_content/summarize.log",
    threads: 1
    resources:
        **gene_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        python {input.script:q} --config {params.config:q} --table {input.table:q} --checkm2 {input.checkm2:q} \
          --genomes {input.genomes:q} --genes {input.genes:q} --kofam {input.kofam:q} --ko-list {input.ko_list:q} \
          --scratch {params.scratch:q} --genome-stats {output.genomes:q} --ko-matrix {output.kos:q} \
          --modules {output.modules:q} --lineages {output.lineages:q} --focal {output.focal:q} > {log:q} 2>&1
        """


rule mycoplasmatales_gene_content:
    input:
        f"{GENE}/summary/genome_stats.tsv",
        f"{GENE}/summary/lineage_modules.tsv",
        f"{GENE}/summary/focal_functions.tsv",
