"""
Phylogenomic placement of the siphonophore Mycoplasmatales MAG species among GTDB
r220 Mycoplasmatales and Mollicutes genomes from cnidarian and planktonic hosts
(config/mycoplasmatales_phylogeny.json; external genomes frozen in
config/mycoplasmatales_external.tsv).

1. Obtain external genomes from NCBI with MD5 checks: fetch_mycoplasmatales_external.
2. Screen external genomes: assess_external_mollicutes_checkm2,
   classify_external_mollicutes_gtdbtk.
3. Select genomes in the order with sufficient completeness: select_phylogeny_genomes.
4. Infer the de novo bac120 tree: infer_mycoplasmatales_tree_gtdbtk.

Next: gene-content comparison of the near-complete genomes uses this placement.
"""

MYCO = f"{WORK}/mycoplasmatales_phylogeny"
with open(ROOT / "config" / "mycoplasmatales_phylogeny.json") as handle:
    MYCO_CONFIG = json.load(handle)
MYCO_EXTERNAL = str(ROOT / MYCO_CONFIG["external_genomes"])


rule fetch_mycoplasmatales_external:
    input:
        script="workflow/scripts/fetch_ncbi_genomes.py",
        table=MYCO_EXTERNAL,
    output:
        genomes=directory(f"{MYCO}/external_genomes"),
        manifest=f"{MYCO}/external_genomes.manifest.tsv",
    log:
        f"{WORK}/logs/mycoplasmatales_phylogeny/fetch_external.log",
    threads: 1
    resources:
        **phase4_resources("download"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --table {input.table:q} --genomes {output.genomes:q} \
          --manifest {output.manifest:q} > {log:q} 2>&1
        """


rule assess_external_mollicutes_checkm2:
    input:
        genomes=f"{MYCO}/external_genomes",
        database=PHASE4_CATALOG["databases"]["checkm2"],
    output:
        f"{MYCO}/external_checkm2.tsv",
    params:
        tools=P4_TOOLS["checkm2"],
        scratch=f"{SCRATCH}/mycoplasmatales_phylogeny/checkm2_tmp",
    log:
        f"{WORK}/logs/mycoplasmatales_phylogeny/checkm2.log",
    threads: 16
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


rule classify_external_mollicutes_gtdbtk:
    input:
        genomes=f"{MYCO}/external_genomes",
        reference=ancient(f"{CATALOG_SCRATCH}/gtdbtk_r220"),
    output:
        bacteria=f"{MYCO}/external_gtdbtk.bac120.summary.tsv",
        archaea=f"{MYCO}/external_gtdbtk.ar53.summary.tsv",
    params:
        tools=P4_TOOLS["gtdbtk"],
        scratch=f"{SCRATCH}/mycoplasmatales_phylogeny/gtdbtk_tmp",
        mash_db=f"{CATALOG_SCRATCH}/gtdbtk_mash/gtdb_r220.msh",
    log:
        f"{WORK}/logs/mycoplasmatales_phylogeny/classify_external.log",
    threads: 16
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
        gtdbtk classify_wf --genome_dir {input.genomes:q} --extension fa --out_dir "$tmp/out" --cpus {threads} \
          --pplacer_cpus 8 --mash_db {params.mash_db:q} >> {log:q} 2>&1
        header=$(head -n 1 "$tmp/out/gtdbtk.bac120.summary.tsv" 2>/dev/null || echo user_genome)
        [ -e "$tmp/out/gtdbtk.ar53.summary.tsv" ] || echo "$header" > "$tmp/out/gtdbtk.ar53.summary.tsv"
        # Top-level summaries are relative symlinks into $tmp; copy their targets.
        cp -L "$tmp/out/gtdbtk.ar53.summary.tsv" {output.archaea:q}.tmp
        cp -L "$tmp/out/gtdbtk.bac120.summary.tsv" {output.bacteria:q}.tmp
        mv {output.archaea:q}.tmp {output.archaea:q}
        mv {output.bacteria:q}.tmp {output.bacteria:q}
        """


rule select_phylogeny_genomes:
    input:
        script="workflow/scripts/select_phylogeny_genomes.py",
        external=MYCO_EXTERNAL,
        genomes=f"{MYCO}/external_genomes",
        checkm2=f"{MYCO}/external_checkm2.tsv",
        gtdbtk=[f"{MYCO}/external_gtdbtk.bac120.summary.tsv", f"{MYCO}/external_gtdbtk.ar53.summary.tsv"],
        mag_catalog=f"{WORK}/phase3_catalog/mags/mag_catalog.tsv",
        mags=f"{WORK}/phase3_catalog/mags/species_representatives",
    output:
        genomes=directory(f"{MYCO}/tree_genomes"),
        labels=f"{MYCO}/tree_labels.tsv",
    params:
        order=MYCO_CONFIG["order"],
        min_completeness=MYCO_CONFIG["min_completeness"],
    log:
        f"{WORK}/logs/mycoplasmatales_phylogeny/select.log",
    threads: 1
    resources:
        **phase4_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --external {input.external:q} --external-genomes {input.genomes:q} \
          --checkm2 {input.checkm2:q} --gtdbtk {input.gtdbtk:q} --mag-catalog {input.mag_catalog:q} \
          --mag-representatives {input.mags:q} --order {params.order} --min-completeness {params.min_completeness} \
          --genomes {output.genomes:q} --labels {output.labels:q} > {log:q} 2>&1
        """


rule infer_mycoplasmatales_tree_gtdbtk:
    input:
        genomes=f"{MYCO}/tree_genomes",
        reference=ancient(f"{CATALOG_SCRATCH}/gtdbtk_r220"),
    output:
        tree=f"{MYCO}/mycoplasmatales.bac120.decorated.tree",
        msa=f"{MYCO}/mycoplasmatales.bac120.msa.fasta.gz",
    params:
        tools=P4_TOOLS["gtdbtk"],
        scratch=f"{SCRATCH}/mycoplasmatales_phylogeny/denovo_tmp",
        taxa_filter=MYCO_CONFIG["gtdbtk_de_novo"]["taxa_filter"],
        outgroup=MYCO_CONFIG["gtdbtk_de_novo"]["outgroup_taxon"],
    log:
        f"{WORK}/logs/mycoplasmatales_phylogeny/de_novo.log",
    benchmark:
        f"{WORK}/benchmarks/mycoplasmatales_phylogeny/de_novo.tsv",
    threads: 16
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
        gtdbtk de_novo_wf --genome_dir {input.genomes:q} --extension fa --bacteria --out_dir "$tmp/out" \
          --taxa_filter {params.taxa_filter} --outgroup_taxon {params.outgroup} --cpus {threads} >> {log:q} 2>&1
        cp -L "$tmp/out/gtdbtk.bac120.decorated.tree" {output.tree:q}.tmp
        cp -L "$tmp/out/align/gtdbtk.bac120.msa.fasta.gz" {output.msa:q}.tmp
        mv {output.tree:q}.tmp {output.tree:q}
        mv {output.msa:q}.tmp {output.msa:q}
        """


rule mycoplasmatales_phylogeny:
    input:
        f"{MYCO}/mycoplasmatales.bac120.decorated.tree",
        f"{MYCO}/tree_labels.tsv",
