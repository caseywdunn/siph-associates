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
        bacterial_data="figures/cohort/bacterial_detections.tsv",
        viral_data="figures/cohort/viral_detections.tsv",
        eukaryote_data="figures/cohort/eukaryote_evidence.tsv",
        flowcell_data="figures/cohort/library_flowcells.tsv",
        incidence_data="figures/cohort/mycoplasmatales_incidence.tsv",
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


rule collection_context_tables:
    input:
        script="scripts/analyze_collection_context.py",
        manifest="manifest.csv",
        bacterial="figures/cohort/bacterial_detections.tsv",
        viral="figures/cohort/viral_detections.tsv",
        eukaryotes="figures/cohort/eukaryote_evidence.tsv",
        flowcells="figures/cohort/library_flowcells.tsv",
        incidence="figures/cohort/mycoplasmatales_incidence.tsv",
    output:
        context="figures/collection_context/library_context.tsv",
        detections="figures/collection_context/associate_detections.tsv",
        host="figures/collection_context/host_detection_summary.tsv",
        depth="figures/collection_context/depth_summary.tsv",
        tissue="figures/collection_context/tissue_summary.tsv",
        dt68="figures/collection_context/dt68_specimens.tsv",
        dt68_summary="figures/collection_context/dt68_depth_summary.tsv",
        regions="figures/collection_context/tissue_region_incidence.tsv",
        host_regions="figures/collection_context/tissue_region_host_incidence.tsv",
        counts="figures/collection_context/summary.json",
        provenance="figures/collection_context/provenance.json",
        source_snapshot="figures/collection_context/provenance/analyze_collection_context.py",
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        partition="day",
    log:
        "logs/figures/collection_context_tables.log",
    shell:
        "/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python {input.script:q} > {log:q} 2>&1"


rule publication_collection_context_figures:
    input:
        script="figures/plot_collection_context.py",
        context="figures/collection_context/library_context.tsv",
        host="figures/collection_context/host_detection_summary.tsv",
        dt68="figures/collection_context/dt68_specimens.tsv",
        tissue="figures/collection_context/tissue_summary.tsv",
        regions="figures/collection_context/tissue_region_incidence.tsv",
        provenance="figures/collection_context/provenance.json",
    output:
        depth="figures/collection_context/depth_host_overview.pdf",
        dt68="figures/collection_context/dt68_depth_context.pdf",
        tissue="figures/collection_context/tissue_context.pdf",
        depth_preview="figures/collection_context/depth_host_overview.png",
        dt68_preview="figures/collection_context/dt68_depth_context.png",
        tissue_preview="figures/collection_context/tissue_context.png",
        provenance="figures/collection_context/plot_manifest.json",
        source_snapshot="figures/collection_context/provenance/plot_collection_context.py",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=5,
        partition="day",
    log:
        "logs/figures/collection_context_plots.log",
    shell:
        "MPLCONFIGDIR=/tmp/siph-associates-matplotlib /gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/python {input.script:q} > {log:q} 2>&1"
