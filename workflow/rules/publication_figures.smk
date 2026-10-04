"""Small publication exports built directly from accepted result tables."""

rule publication_cohort_figures:
    input:
        script="figures/build_cohort.py",
        manifest="manifest.csv",
        assembly="config/phase3_assembly.tsv",
        bacterial=f"{P6V}/grades/bacterial_grades.tsv",
        viral=f"{P6V}/grades/viral_grades.tsv",
        flowcells=f"{P6V}/primary/library_flowcells.tsv",
        eukaryotes=f"{WORK}/eukaryote_gate/eukaryote_grades.tsv",
        mags=f"{WORK}/phase3_catalog/mags/mag_catalog.tsv",
        catalog=f"{P6_CATALOG}/bacterial_catalog.manifest.tsv",
    output:
        sampling="figures/cohort/sampling.pdf",
        incidence="figures/cohort/mycoplasmatales_incidence.pdf",
        eukaryotes="figures/cohort/eukaryote_evidence.pdf",
        counts="figures/cohort/summary.json",
        provenance="figures/cohort/source_checksums.tsv",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=15,
        partition="day",
    log:
        "logs/figures/cohort.log",
    shell:
        "/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python {input.script:q} > {log:q} 2>&1"


rule publication_mycoplasmatales_figures:
    input:
        script="workflow/scripts/plot_mycoplasmatales_figures.py",
        focal=f"{WORK}/mycoplasmatales_gene_content/summary/focal_functions.tsv",
        genomes=f"{WORK}/mycoplasmatales_gene_content/summary/genome_stats.tsv",
        tree=f"{WORK}/mycoplasmatales_phylogeny/mycoplasmatales.bac120.decorated.tree",
        labels=f"{WORK}/mycoplasmatales_phylogeny/tree_labels.tsv",
    output:
        functions="figures/mycoplasmatales/gene_functions.pdf",
        tree="figures/mycoplasmatales/phylogeny_focus.pdf",
        provenance="figures/mycoplasmatales/manifest.json",
    threads: 1
    resources:
        mem_mb=4000,
        runtime=15,
        partition="day",
    log:
        "logs/figures/mycoplasmatales.log",
    shell:
        "/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python {input.script:q} > {log:q} 2>&1"
