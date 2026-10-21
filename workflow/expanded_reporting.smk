"""Downstream reporting for the completed 213-library v3 mappings.

Imports the shared evidence/statistics rules. Advertised targets reuse completed
assemblies, catalogs and cohort mappings; whole-library read controls are disabled.
Historical v2 read-handling controls retain their original scope.
"""
import csv
import json
from pathlib import Path

configfile: "config/expanded_reporting.json"
ROOT = Path.cwd()
WORK, SCRATCH = config['work_root'], config['scratch_root']
PY = config['python']
with open(config['samples']) as handle:
    sample_rows = list(csv.DictReader(handle, delimiter='\t'))
SAMPLES = {r['sample_id']: r for r in sample_rows}
SAMPLE_IDS = sorted(SAMPLES)
ASSEMBLY_IDS = sorted(s for s, r in SAMPLES.items() if r['assembly_eligible'] == 'true')
if len(SAMPLES) != len(sample_rows) or len(SAMPLES) != config['expected_samples']:
    raise WorkflowError('Expanded reporting requires exactly 213 unique libraries')
PHASE5 = json.loads(Path('config/phase5_mapping.json').read_text())
PHASE5['catalog_version'] = config['catalog_version']
PHASE6_CONFIG_PATH = ROOT / 'config/phase6_analysis.json'
PHASE6 = json.loads(PHASE6_CONFIG_PATH.read_text())
PHASE6_SENSITIVITY = []

def phase5_resources(name):
    resources = {k: v for k, v in PHASE5['resources'][name].items() if k != 'threads'}
    resources['partition'] = 'ycga'
    return resources

def disabled_reads(*args):
    raise WorkflowError('Whole-library read controls are outside this reporting entry point')
raw_r1 = raw_r2 = phase6_sensitivity_reads = phase6_index = disabled_reads

wildcard_constraints:
    sample='[A-Za-z0-9_.-]+',

rule reporting_metadata:
    input:
        code=['scripts/reporting_manifest.py', 'scripts/export_sample_metadata.py'],
        sources=['manifest.csv', 'config/samples.tsv', 'config/publication_host_names.tsv',
                 'data/metadata/sample_metadata_updates.tsv', 'data/metadata/publication_metadata_updates.tsv',
                 'manko2026_manifest.csv', 'data/metadata/manko2026/sample_metadata.tsv',
                 'config/post_expansion_manifest.csv', config['samples']],
    output:
        manifest=config['manifest'],
        tables=expand(f'{WORK}/metadata/{{name}}', name=['sample_metadata.tsv',
            'sample_metadata_provenance.tsv', 'sample_metadata_missing.tsv', 'sample_metadata_export.json']),
    params: directory=f'{WORK}/metadata'
    resources: mem_mb=2000, runtime=10, partition='ycga'
    log: f'{WORK}/logs/metadata.log'
    shell: '{PY:q} scripts/reporting_manifest.py --output-dir {params.directory:q} > {log:q} 2>&1'

include: 'rules/phase6_analysis.smk'

rule expanded_statistics:
    input: rules.phase6_analysis.input
    default_target: True

EUKARYOTE_CONFIG_PATH = ROOT / 'config/eukaryote_gate.json'
EUKARYOTE = json.loads(EUKARYOTE_CONFIG_PATH.read_text())
include: 'rules/eukaryote_gate.smk'

rule expanded_eukaryotes:
    input: rules.eukaryote_gate.input

rule expanded_eukaryote_evidence:
    input:
        f'{EUK}/read_lineages.tsv',
        f'{EUK}/assembled_nonhost.tsv',
        f'{EUK}/assembled_nonhost.fasta',

include: 'rules/publication_figures.smk'

rule expanded_screen_qc:
    input:
        code='scripts/expanded_screening_qc.py',
        config='config/reporting_sources.json',
        baseline='data/results/aggregation/cohort/library_qc.tsv',
        metadata='data/metadata/manko2026/sample_metadata.tsv',
        completed='data/results/bargmannia/2026-10-19/stages/phase3_assembly.done',
        reports=expand('data/results/bargmannia/2026-10-19/screens/bracken/{sample}.{rank}.tsv',
                       sample=[s for s in SAMPLE_IDS if s.startswith('Manko2026__')], rank=['G', 'S']),
    output: PUBLICATION['screen_qc']
    resources: mem_mb=2000, runtime=10, partition='ycga'
    log: f'{WORK}/logs/screens/qc.log'
    shell: '{PY:q} {input.code:q} > {log:q} 2>&1'

rule expanded_nomination_evidence:
    input:
        code='scripts/phase4_nomination_evidence.py',
        config='config/reporting_sources.json',
        nominations=PUBLICATION['aggregation_root'] + '/candidate_nominations.tsv',
        mags=PUBLICATION['mag_catalog'],
        taxonomy='/gpfs/ycga/work/dunn/cwd7/databases/gtdbtk_r220/release220/taxonomy/gtdb_taxonomy.tsv',
    output:
        expand(PUBLICATION['nomination_evidence'] + '/{artifact}.tsv',
               artifact=['sylph_nominations', 'gtdb_corroboration', 'gtdb_uncorroborated', 'contaminant_genera', 'mag_overlap'])
    resources: mem_mb=4000, runtime=15, partition='ycga'
    log: f'{WORK}/logs/screens/nominations.log'
    shell: '{PY:q} {input.code:q} --reporting-config {input.config:q} > {log:q} 2>&1'

rule expanded_initial_ssu_evidence:
    input:
        code='scripts/extract_eukaryote_read_evidence.py',
        baseline='data/results/eukaryote_gate/read_evidence.tsv',
        archives=expand('data/results/bargmannia/2026-10-19/screens/phyloflash/{sample}.tar.gz',
                        sample=[s for s in SAMPLE_IDS if s.startswith('Manko2026__')]),
    output: f'{EUK}/read_evidence.tsv'
    resources: mem_mb=16000, runtime=120, partition='ycga'
    log: f'{WORK}/logs/screens/initial_ssu.log'
    shell: '{PY:q} {input.code:q} --source data/results/bargmannia/2026-10-19/screens/phyloflash --baseline {input.baseline:q} --output {output:q} > {log:q} 2>&1'

rule expanded_ssu_controls:
    input:
        code='scripts/calibrate_eukaryote_reads.py',
        units='workflow/scripts/eukaryote_units.py',
        reads=f'{EUK}/read_lineages.tsv',
        classifications=[eukaryote_ssu(s, 'ssu_classification.csv') for s in SAMPLE_IDS],
    output:
        expand(f'{EUK}/calibration/{{name}}.tsv', name=['negative_controls',
            'assembled_positives', 'threshold_grid', 'detections_above_background_ge97',
            'detections_above_background_ge99'])
    resources: mem_mb=2000, runtime=10, partition='ycga'
    log: f'{WORK}/logs/screens/ssu_controls.log'
    shell: '{PY:q} {input.code:q} --gate-root {EUK:q} --ssu-dir data/results/eukaryote_gate/ssu --ssu-dir {EUK:q}/ssu > {log:q} 2>&1'

rule expanded_screening_summary:
    input:
        code='scripts/summarize_screening_results.py',
        config='config/reporting_sources.json',
        policies=['config/eukaryote_gate.json', 'config/phase6_analysis.json'],
        catalog=PUBLICATION['catalog_manifest'],
        raw_nominations=PUBLICATION['aggregation_root'] + '/candidate_nominations.tsv',
        read_lineages=f'{EUK}/read_lineages.tsv',
        qc=rules.expanded_screen_qc.output,
        nominations=rules.expanded_nomination_evidence.output,
        initial=rules.expanded_initial_ssu_evidence.output,
        controls=rules.expanded_ssu_controls.output,
        grades=f'{P6V}/grades/bacterial_grades.tsv',
        eukaryotes=f'{EUK}/eukaryote_grades.tsv',
        sensitivity=f'{P6V}/sensitivity/leave_one_study_out.tsv',
    output:
        expand('figures/screening/{name}.tsv', name=['method_summary',
            'sylph_nomination_grades', 'ssu_control_comparison', 'sylph_sources',
            'genus_corroboration', 'bacterial_detection_overlap',
            'mycoplasmatales_screen_support', 'eukaryote_detection_overlap']),
        'figures/screening/summary.json', 'figures/screening/values.tex',
        'figures/screening/provenance.json',
    resources: mem_mb=8000, runtime=20, partition='ycga'
    log: f'{WORK}/logs/screens/summary.log'
    shell: '{PY:q} {input.code:q} > {log:q} 2>&1'

rule expanded_screen_evidence:
    input:
        rules.expanded_screen_qc.output,
        rules.expanded_nomination_evidence.output,
        rules.expanded_initial_ssu_evidence.output,

rule expanded_overview:
    input:
        script='figures/plot_associate_overview.py',
        settings='config/reporting_sources.json',
        cohort=rules.publication_cohort_figures.output,
        context=rules.collection_context_tables.output,
    output:
        expand('figures/overview/{name}.{extension}',
               name=['study_overview', 'microbial_host_distribution', 'eukaryote_evidence'],
               extension=['pdf', 'png']),
        expand('figures/overview/{name}.tsv', name=['study_sampling',
            'microbial_host_distribution', 'microbial_display_assignments',
            'eukaryote_display', 'domain_summary']),
        'figures/overview/manifest.json',
    resources: mem_mb=4000, runtime=15, partition='ycga'
    log: f'{WORK}/logs/overview.log'
    shell: 'MPLCONFIGDIR=/tmp/siph-associates-matplotlib {PY:q} {input.script:q} > {log:q} 2>&1'

rule expanded_supporting_tables:
    input:
        scripts=['scripts/export_domain_review_statistics.py', 'scripts/summarize_dalila_revision.py', 'scripts/summarize_mdb_revision.py'],
        settings='config/reporting_sources.json',
        primary=rules.validate_phase6.output,
        statistics=rules.validate_phase6.input,
        genes=expand(PUBLICATION['gene_content_root'] + '/summary/{table}.tsv',
                     table=['genome_stats', 'focal_functions', 'module_completeness', 'module_completeness_strict']),
        assemblies=PUBLICATION['assembly_summary'],
        mags=PUBLICATION['mag_catalog'],
        sources=config['samples'],
        cohort=rules.publication_cohort_figures.output,
        context=rules.collection_context_tables.output,
    output:
        expand('figures/catalog_review/{name}', name=['statistics.json',
            'statistics.tex', 'statistics_provenance.json', 'threshold_rows.tex',
            'threshold_table.tex']),
        expand('figures/dalila_revision/{name}', name=['assembly_routes.tsv',
            'clade_a_regions.tsv', 'host_species_summary.tsv', 'recorded_life_stages.tsv',
            'summary.json', 'values.tex', 'provenance.json']),
        expand('figures/mdb_revision/{name}', name=['amino_acid_modules.tsv',
            'cooccurrence.tsv', 'cooccurrence_results.tex', 'cooccurrence_strata.tsv',
            'cooccurrence_table.tex', 'defense_chitin_genes.tsv', 'genome_quality.tsv',
            'nutrition_defense_table.tex', 'phage_summary.tex', 'specimens.tsv', 'provenance.json']),
    resources: mem_mb=4000, runtime=20, partition='ycga'
    log: f'{WORK}/logs/supporting_tables.log'
    shell: 'for script in {input.scripts:q}; do {PY:q} "$script"; done > {log:q} 2>&1'

rule expanded_publication:
    input:
        rules.expanded_statistics.input,
        rules.expanded_eukaryotes.input,
        rules.publication_collection_context_figures.output,
        rules.expanded_overview.output,
        rules.expanded_supporting_tables.output,
        rules.expanded_screening_summary.output,
        'figures/fish_12s/summary.json',
        'figures/expanded_review/summary.json',
        'figures/host_depletion/mapping_stats.tsv',

rule expanded_fish_summary:
    input:
        code='scripts/fish_12s/finalize_report.py',
        helpers=['scripts/fish_12s/pilot.py', 'scripts/fish_12s/summarize_cohort.py'],
        settings=['config/fish_12s_full_reads.tsv', 'config/fish_12s_host_depleted_cohort.json',
                  'data/metadata/manko2026/sample_metadata.tsv'],
        metadata=rules.reporting_metadata.output.tables,
        existing=expand('data/results/fish_12s/final_review/2026-10-19/host_depleted/{name}',
                        name=['run.json', 'sample_summary.tsv', 'pair_review.tsv']),
        originals=expand('data/results/fish_12s/full_reads/2026-10-18/review/{sample}/{name}',
                         sample=[r['analysis_sample_id'] for r in csv.DictReader(open('config/fish_12s_full_reads.tsv'), delimiter='\t')],
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

rule compare_expanded_publication:
    input:
        script='scripts/compare_expanded_reporting.py',
        units='workflow/scripts/eukaryote_units.py',
        policy='config/eukaryote_gate.json',
        settings='config/reporting_sources.json',
        lineages=f'{EUK}/read_lineages.tsv',
        specimen_metadata='data/metadata/manko2026/sample_metadata.tsv',
        cohort=rules.publication_cohort_figures.output,
        baseline=f'{WORK}/baseline/figures-44ade17.tar.gz',
        catalogs=['figures/catalog_qc/catalog_v3/bacterial_catalog.manifest.tsv',
                  'figures/catalog_qc/catalog_v3/mag_catalog.tsv',
                  'figures/catalog_qc/catalog_v3/catalog_id_crosswalk.tsv'],
        samples=config['samples'],
    output:
        'figures/expanded_review/summary.json',
        'figures/expanded_review/detection_changes.tsv',
        'figures/expanded_review/original_metadata_changes.tsv',
        'figures/expanded_review/retained_focal_genomes.tsv',
        'figures/expanded_review/ssu_background_flags.tsv',
        'figures/expanded_review/background_sensitivity.tsv',
    resources: mem_mb=4000, runtime=10, partition='ycga'
    log: f'{WORK}/logs/comparison.log'
    shell: '{PY:q} {input.script:q} > {log:q} 2>&1'

config['host_depletion_require_complete'] = True
include: 'rules/host_depletion_report.smk'
