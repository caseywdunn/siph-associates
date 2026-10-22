"""Reporting-only interpretation and metadata review of completed cohort evidence.

No rule in this entry point reads FASTQs, assembles sequences or maps reads.
"""
import json
from pathlib import Path

configfile: "config/coauthor_reporting.json"
WORK = config["work_root"]
P6V = f"{WORK}/phase6_analysis/v3"
PY = config["python"]
BASE = config["baseline_root"]
CAT = config["catalog_output"]
META = f"{WORK}/metadata/manifest.csv"
REUSED_STATS = [f"{suffix}/{name}.tsv" for suffix in
    ["primary", "sensitivity/breadth_5pct", "sensitivity/breadth_20pct"]
    for name in ["physalia_design", "physalia_models", "physalia_permanova", "physalia_region_permutation"]]
REUSED_SENSITIVITY = [f"sensitivity/{name}.tsv" for name in ["threshold_sensitivity", "leave_one_study_out"]]

rule coauthor_reporting_metadata:
    input:
        code=["scripts/reporting_manifest.py", "scripts/export_sample_metadata.py", "workflow/scripts/phase6_metadata.py"],
        metadata=["manifest.csv", "config/samples.tsv", "config/publication_host_names.tsv",
                  "data/metadata/sample_metadata_updates.tsv", "data/metadata/publication_metadata_updates.tsv",
                  "manko2026_manifest.csv", "data/metadata/manko2026/sample_metadata.tsv",
                  "config/post_expansion_manifest.csv", config["samples"]],
    output:
        manifest=META,
        flow=f"{P6V}/primary/library_flowcells.tsv",
    params: directory=f"{WORK}/metadata"
    log: f"{WORK}/logs/metadata.log"
    resources: mem_mb=2000, runtime=10
    shell:
        """
        {PY:q} scripts/reporting_manifest.py --output-dir {params.directory:q} > {log:q} 2>&1
        {PY:q} workflow/scripts/phase6_metadata.py --manifest {output.manifest:q} --flowcell-policy single_flowcell --output {output.flow:q} >> {log:q} 2>&1
        """

rule review_viral_classification:
    input:
        script="scripts/prepare_coauthor_reporting.py",
        configuration="config/coauthor_reporting.json",
        changes=config["classifications"],
        completed=config["baseline_validation"],
        oldflow=f"{BASE}/primary/library_flowcells.tsv",
        newflow=rules.coauthor_reporting_metadata.output.flow,
        grades=expand(f"{BASE}/grades/{{kind}}_grades.tsv", kind=["bacterial", "viral"]),
        statistics=[f"{BASE}/{path}" for path in REUSED_STATS + REUSED_SENSITIVITY],
        manifests=expand(config["catalog_source"] + "/{name}", name=["bacterial_catalog.manifest.tsv", "mag_catalog.tsv", "host_links.tsv", "catalog_id_crosswalk.tsv", "viral_associate.manifest.tsv", "viral_endogenous_candidate.manifest.tsv"]),
    output:
        prepared=f"{WORK}/prepared.json",
        bacterial=f"{P6V}/grades/bacterial_grades.tsv",
        viral=f"{P6V}/grades/viral_grades.tsv",
        statistics=[f"{P6V}/{path}" for path in REUSED_STATS + REUSED_SENSITIVITY],
        catalog=expand(CAT + "/{name}", name=["bacterial_catalog.manifest.tsv", "mag_catalog.tsv", "host_links.tsv", "catalog_id_crosswalk.tsv", "viral_associate.manifest.tsv", "viral_endogenous_candidate.manifest.tsv", "reclassified_detections.tsv", "summary.json", "provenance.json"]),
    log: f"{WORK}/logs/classification.log"
    resources: mem_mb=4000, runtime=10
    shell: "{PY:q} {input.script:q} --config {input.configuration:q} > {log:q} 2>&1"

rule reviewed_flowcell_tests:
    input:
        script="workflow/scripts/test_contamination_flowcell.py",
        analysis="config/phase6_analysis.json",
        bacterial=rules.review_viral_classification.output.bacterial,
        flow=rules.coauthor_reporting_metadata.output.flow,
    output:
        cohort=f"{P6V}/review_tests/{{threshold}}/contamination_tests.tsv",
        nanomia=f"{P6V}/review_tests/{{threshold}}/nanomia_contamination_tests.tsv",
    params: grade=lambda w: "grade" if w.threshold == "10" else f"grade_at_{w.threshold}pct"
    wildcard_constraints: threshold="5|10|20"
    log: f"{WORK}/logs/flowcells_{{threshold}}.log"
    resources: mem_mb=4000, runtime=15
    shell:
        """
        {PY:q} {input.script:q} --analysis {input.analysis:q} --grades {input.bacterial:q} --metadata {input.flow:q} --grade-column {params.grade:q} --output {output.cohort:q} > {log:q} 2>&1
        {PY:q} {input.script:q} --analysis {input.analysis:q} --grades {input.bacterial:q} --metadata {input.flow:q} --grade-column {params.grade:q} --subset nanomia --output {output.nanomia:q} >> {log:q} 2>&1
        """

rule reviewed_incidence:
    input:
        script="workflow/scripts/summarize_incidence.py",
        bacterial=rules.review_viral_classification.output.bacterial,
        viral=rules.review_viral_classification.output.viral,
        contamination=f"{P6V}/review_tests/{{threshold}}/contamination_tests.tsv",
        manifest=META,
        samples=config["samples"],
        links=f"{CAT}/host_links.tsv",
    output:
        bacterial=f"{P6V}/review_tests/{{threshold}}/bacterial_incidence.tsv",
        viral=f"{P6V}/review_tests/{{threshold}}/viral_incidence.tsv",
        grades=f"{P6V}/review_tests/{{threshold}}/grade_summary.tsv",
        links=f"{P6V}/review_tests/{{threshold}}/phage_host_links.tsv",
    params: grade=lambda w: "grade" if w.threshold == "10" else f"grade_at_{w.threshold}pct"
    wildcard_constraints: threshold="5|10|20"
    log: f"{WORK}/logs/incidence_{{threshold}}.log"
    resources: mem_mb=4000, runtime=10
    shell:
        """
        {PY:q} {input.script:q} --bacterial {input.bacterial:q} --viral {input.viral:q} --contamination {input.contamination:q} --manifest {input.manifest:q} --samples {input.samples:q} --links {input.links:q} --grade-column {params.grade:q} --out-bacterial {output.bacterial:q} --out-viral {output.viral:q} --out-grades {output.grades:q} --out-links {output.links:q} > {log:q} 2>&1
        """

rule reviewed_reporting_validation:
    input:
        code="scripts/validate_coauthor_reporting.py",
        configuration="config/coauthor_reporting.json",
        prepared=rules.review_viral_classification.output,
        tables=expand(f"{P6V}/review_tests/{{threshold}}/{{name}}.tsv", threshold=[5,10,20], name=["contamination_tests", "nanomia_contamination_tests", "bacterial_incidence", "viral_incidence", "grade_summary", "phage_host_links"]),
        metadata=rules.coauthor_reporting_metadata.output,
    output:
        done=f"{WORK}/validated.done",
        primary=expand(f"{P6V}/primary/{{name}}.tsv", name=["nanomia_contamination_tests", "bacterial_incidence", "viral_incidence", "grade_summary", "phage_host_links"]),
        contamination=f"{P6V}/grades/contamination_tests.tsv",
        thresholds=expand(f"{P6V}/sensitivity/breadth_{{threshold}}pct/{{name}}.tsv", threshold=[5,20], name=["contamination_tests", "nanomia_contamination_tests", "bacterial_incidence", "viral_incidence", "grade_summary", "phage_host_links"]),
    log: f"{WORK}/logs/validation.log"
    resources: mem_mb=4000, runtime=10
    shell: "{PY:q} {input.code:q} --config {input.configuration:q} > {log:q} 2>&1"

rule reviewed_publication_tables:
    input:
        validated=rules.reviewed_reporting_validation.output,
        scripts=["figures/build_cohort.py", "scripts/analyze_collection_context.py", "figures/plot_collection_context.py", "figures/plot_associate_overview.py", "scripts/export_domain_review_statistics.py", "scripts/summarize_dalila_revision.py", "scripts/summarize_mdb_revision.py", "scripts/compare_expanded_reporting.py", "scripts/summarize_screening_results.py", "scripts/fish_12s/finalize_report.py"],
        settings="config/reporting_sources.json",
    output:
        done="figures/coauthor_review/publication.done",
        cohort=expand("figures/cohort/{name}", name=["sample_metadata.tsv", "bacterial_detections.tsv", "viral_detections.tsv", "eukaryote_evidence.tsv", "library_flowcells.tsv", "mycoplasmatales_incidence.tsv", "summary.json"]),
        summaries=["figures/dalila_revision/summary.json", "figures/catalog_review/statistics.json", "figures/expanded_review/summary.json", "figures/collection_context/summary.json", "figures/fish_12s/summary.json", "figures/screening/summary.json"],
    log: f"{WORK}/logs/publication.log"
    resources: mem_mb=8000, runtime=20
    shell:
        """
        export MPLCONFIGDIR=/tmp/siph-associates-matplotlib
        for script in {input.scripts:q}; do {PY:q} "$script"; done > {log:q} 2>&1
        touch {output.done:q}
        """

rule coauthor_publication:
    input:
        rules.reviewed_reporting_validation.output,
        rules.reviewed_publication_tables.output,
    default_target: True
