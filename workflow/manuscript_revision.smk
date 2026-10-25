"""Bounded revision of completed mapping evidence; no read mapping is scheduled.

Assembly support is realigned before regrading and contamination-aware reporting.
The accepted per-library read coverage and assembly files remain frozen inputs.
"""
import csv
import json
from pathlib import Path

configfile: "config/manuscript_revision.json"
ROOT = Path.cwd()
WORK, SCRATCH = config["work_root"], config["scratch_root"]
PY = config["python"]
with open(config["samples"]) as handle:
    sample_rows = list(csv.DictReader(handle, delimiter="\t"))
SAMPLES = {row["sample_id"]: row for row in sample_rows}
SAMPLE_IDS = sorted(SAMPLES)
ASSEMBLY_IDS = sorted(s for s, row in SAMPLES.items() if row["assembly_eligible"] == "true")
assert len(SAMPLES) == 213 and len(ASSEMBLY_IDS) == 213
PHASE5 = json.loads(Path("config/phase5_mapping.json").read_text())
PHASE5["catalog_version"] = "v3"
PHASE6_CONFIG_PATH = ROOT / "config/phase6_analysis.json"
PHASE6 = json.loads(PHASE6_CONFIG_PATH.read_text())
PHASE6_SENSITIVITY = []

def phase5_resources(name):
    resources = {key: value for key, value in PHASE5["resources"][name].items() if key != "threads"}
    resources["partition"] = "ycga"
    return resources

def disabled_reads(*args):
    raise WorkflowError("Read remapping is outside the authorized revision")
raw_r1 = raw_r2 = phase6_sensitivity_reads = phase6_index = disabled_reads

include: "rules/phase6_analysis.smk"

rule revision_assembly_support:
    input: expand(f"{P6V}/assembly_support/{{sample}}.paf", sample=ASSEMBLY_IDS)

rule revision_statistics:
    input: rules.phase6_analysis.input
    default_target: True


rule revision_catalog:
    input:
        script="scripts/prepare_revision_catalog.py",
        changes="config/viral_catalog_classifications.tsv",
        sources=expand("figures/coauthor_review/catalog/{name}", name=["bacterial_catalog.manifest.tsv", "mag_catalog.tsv", "catalog_id_crosswalk.tsv"]),
        viruses=expand("figures/catalog_qc/catalog_v3/viral_{kind}.manifest.tsv", kind=["associate", "endogenous_candidate"]),
    output:
        expand("figures/manuscript_revision/catalog/{name}", name=["bacterial_catalog.manifest.tsv", "mag_catalog.tsv", "catalog_id_crosswalk.tsv", "viral_associate.manifest.tsv", "viral_endogenous_candidate.manifest.tsv", "summary.json", "provenance.json"])
    resources: mem_mb=2000, runtime=5, partition="ycga"
    shell: "{PY:q} {input.script:q}"

include: "rules/publication_figures.smk"

rule revision_cohort:
    input: rules.publication_cohort_figures.output

rule revision_descriptive_reporting:
    input:
        cohort=rules.publication_cohort_figures.output,
        context=rules.publication_collection_context_figures.output,
        statistics=rules.phase6_analysis.input,
        scripts=["figures/plot_associate_overview.py", "scripts/export_domain_review_statistics.py", "scripts/summarize_dalila_revision.py", "scripts/summarize_screening_results.py", "scripts/summarize_contamination_revision.py", "scripts/summarize_human_cutibacterium.py"],
        settings="config/reporting_sources.json",
        decisions=["config/excluded_libraries.tsv", "config/contaminant_taxa.tsv"],
    output:
        "figures/overview/manifest.json",
        "figures/catalog_review/statistics.json",
        "figures/dalila_revision/summary.json",
        "figures/screening/summary.json",
        "figures/contamination/provenance.json",
        "figures/human_cutibacterium/provenance.json",
    resources: mem_mb=8000, runtime=30, partition="ycga"
    log: f"{WORK}/logs/descriptive_reporting.log"
    shell:
        "for script in {input.scripts:q}; do {PY:q} \"$script\"; done > {log:q} 2>&1"

# Refresh reports from completed fish screens; never rescan reads.
FISH_HALISTEMMA_CONFIG = 'config/fish_12s_halistemma_host_depleted.json'
FISH_HALISTEMMA = json.loads(Path(FISH_HALISTEMMA_CONFIG).read_text())
with open(FISH_HALISTEMMA['manifest']) as handle:
    FISH_HALISTEMMA_IDS = {r['analysis_sample_id'] for r in csv.DictReader(handle, delimiter='\t')}

rule revision_fish_summary:
    input:
        code='scripts/fish_12s/finalize_report.py',
        helpers=['scripts/fish_12s/pilot.py', 'scripts/fish_12s/summarize_cohort.py'],
        settings=['config/fish_12s_full_reads.tsv', 'config/fish_12s_host_depleted_cohort.json',
                  FISH_HALISTEMMA_CONFIG, FISH_HALISTEMMA['manifest'],
                  'data/metadata/manko2026/sample_metadata.tsv'],
        metadata="figures/cohort/sample_metadata.tsv",
        existing=expand('data/results/fish_12s/final_review/2026-10-19/host_depleted/{name}',
                        name=['run.json', 'sample_summary.tsv', 'pair_review.tsv']),
        originals=expand('data/results/fish_12s/full_reads/2026-10-18/review/{sample}/{name}',
                         sample=[r['analysis_sample_id'] for r in csv.DictReader(open('config/fish_12s_full_reads.tsv'), delimiter='\t')
                                 if r['analysis_sample_id'] not in FISH_HALISTEMMA_IDS],
                         name=['run.json', 'sample_summary.tsv', 'pair_review.tsv']),
        replacement=expand(FISH_HALISTEMMA['run'] + '/review/{name}',
                           name=['run.json', 'sample_summary.tsv', 'pair_review.tsv']),
        additions=expand('data/results/bargmannia/2026-10-19/fish_12s/review/{sample}/{name}',
                         sample=[s for s in SAMPLE_IDS if s.startswith('Manko2026__')],
                         name=['run.json', 'sample_summary.tsv', 'pair_review.tsv']),
    output:
        'figures/fish_12s/summary.json', 'figures/fish_12s/sample_summary.tsv',
        'figures/fish_12s/positive_pairs.tsv', 'figures/fish_12s/provenance.json',
        'figures/fish_12s/values.tex',
    resources: mem_mb=4000, runtime=15, partition='ycga'
    log: f'{WORK}/logs/fish_summary.log'
    shell: '{PY:q} {input.code:q} > {log:q} 2>&1'


include: "rules/coauthor_metadata.smk"

rule revision_validation:
    input:
        script="scripts/validate_manuscript_revision.py",
        descriptive=rules.revision_descriptive_reporting.output,
        targeted=rules.record_coauthor_metadata.output,
    output:
        "figures/manuscript_revision/validation.json",
        "figures/manuscript_revision/focal_test_comparison.tsv",
        "figures/manuscript_revision/community_profile_changes.tsv",
    resources: mem_mb=4000, runtime=10, partition="ycga"
    shell: "{PY:q} {input.script:q}"

rule revision_publication:
    input: rules.revision_validation.output
