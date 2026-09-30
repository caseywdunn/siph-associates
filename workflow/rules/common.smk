import json


def split_paths(value):
    return [path for path in value.split(";") if path]


def source_r1(wildcards):
    return split_paths(SAMPLES[wildcards.sample]["r1_paths"])


def source_r2(wildcards):
    return split_paths(SAMPLES[wildcards.sample]["r2_paths"])


def raw_r1(wildcards):
    if config["fixture_enabled"]:
        return [f"{FIXTURE_ROOT}/raw/{wildcards.sample}_R1.fastq.gz"]
    return source_r1(wildcards)


def raw_r2(wildcards):
    if config["fixture_enabled"]:
        return [f"{FIXTURE_ROOT}/raw/{wildcards.sample}_R2.fastq.gz"]
    return source_r2(wildcards)


def host_route(wildcards):
    return SAMPLES[wildcards.sample]["host_route"]


def host_index_inputs(wildcards):
    route = host_route(wildcards)
    if route == "none":
        return []
    prefix = f"{WORK}/reference_indices/{route}/host"
    return [prefix + suffix for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")]


def reference_source(wildcards):
    if config["fixture_enabled"]:
        return f"{FIXTURE_ROOT}/references/{wildcards.route}.fasta"
    return config["databases"]["host_references"][wildcards.route]


PHASE3_CONFIG_PATH = ROOT / "config" / "phase3_benchmark.json"
with PHASE3_CONFIG_PATH.open() as handle:
    PHASE3 = json.load(handle)
PHASE3_PANEL_PATH = ROOT / PHASE3["panel"]
with PHASE3_PANEL_PATH.open(newline="") as handle:
    PHASE3_ROWS = list(csv.DictReader(handle, delimiter="\t")) if config.get("mode") == "production" else []
PHASE3_BY_ID = {row["benchmark_id"]: row for row in PHASE3_ROWS}
PHASE3_IDS = sorted(PHASE3_BY_ID)
PHASE3_SAMPLES = sorted({row["sample_id"] for row in PHASE3_ROWS})

if len(PHASE3_BY_ID) != len(PHASE3_ROWS):
    raise WorkflowError("Phase-3 benchmark IDs are not unique")
if config.get("mode") == "production" and any(row["sample_id"] not in SAMPLES for row in PHASE3_ROWS):
    raise WorkflowError("Phase-3 benchmark contains a sample absent from the production table")


def phase3_row(wildcards):
    try:
        return PHASE3_BY_ID[wildcards.benchmark_id]
    except KeyError:
        raise WorkflowError(f"unknown Phase-3 benchmark ID: {wildcards.benchmark_id}")


def phase3_sample(wildcards):
    return phase3_row(wildcards)["sample_id"]


def phase3_strategy(wildcards):
    return phase3_row(wildcards)["strategy"]


def phase3_read_input(wildcards, mate):
    row = phase3_row(wildcards)
    sample, strategy = row["sample_id"], row["strategy"]
    if strategy == "kraken_nominated":
        return f"{WORK}/screens/kraken/{sample}.classified_{mate}.fastq.gz"
    directory = "fixed" if strategy == "fixed_effort_trimmed" else "reference"
    return f"{SCRATCH}/phase3_benchmark/inputs/{directory}/{wildcards.benchmark_id}_R{mate}.fastq.gz"


def phase3_read_provenance(wildcards):
    row = phase3_row(wildcards)
    if row["strategy"] == "kraken_nominated":
        return f"{WORK}/provenance/kraken_bracken/{row['sample_id']}.json"
    directory = "fixed" if row["strategy"] == "fixed_effort_trimmed" else "reference"
    return f"{WORK}/provenance/phase3_benchmark/inputs/{directory}/{wildcards.benchmark_id}.json"


def phase3_reference_inputs(wildcards):
    route = SAMPLES[phase3_sample(wildcards)]["host_route"]
    if route == "none":
        return []
    prefix = f"{WORK}/reference_indices/{route}/host"
    return [prefix + suffix for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")]


def phase3_metric_paths(kind):
    return expand(f"{WORK}/phase3_benchmark/{kind}/{{benchmark_id}}.json", benchmark_id=PHASE3_IDS)


PHASE3_COHORT_CONFIG_PATH = ROOT / "config" / "phase3_cohort.json"
with PHASE3_COHORT_CONFIG_PATH.open() as handle:
    PHASE3_COHORT = json.load(handle)


def phase3_cohort_strategy(wildcards):
    return PHASE3_COHORT["strategy_by_host_route"][host_route(wildcards)]


def phase3_cohort_input_path(wildcards, mate):
    return f"{WORK}/phase3_cohort/inputs/{wildcards.sample}_R{mate}.fastq.gz"


def phase3_cohort_input_metrics(wildcards):
    return f"{WORK}/phase3_cohort/inputs/{wildcards.sample}.json"


PHASE3_ASSEMBLY_CONFIG_PATH = ROOT / "config" / "phase3_assembly.json"
with PHASE3_ASSEMBLY_CONFIG_PATH.open() as handle:
    PHASE3_ASSEMBLY = json.load(handle)
PHASE3_ASSEMBLY_MEMBERSHIP_PATH = ROOT / PHASE3_ASSEMBLY[
    "fixture_membership" if config["fixture_enabled"] else "membership"
]
with PHASE3_ASSEMBLY_MEMBERSHIP_PATH.open(newline="") as handle:
    PHASE3_ASSEMBLY_ROWS = {row["sample_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
if set(PHASE3_ASSEMBLY_ROWS) != set(SAMPLE_IDS):
    raise WorkflowError("Phase-3 assembly membership does not exactly match the sample table")
ASSEMBLY_IDS = sorted(s for s, row in PHASE3_ASSEMBLY_ROWS.items() if row["assembly_eligible"] == "true")


def phase3_assembly_resource(key):
    def resource(wildcards):
        pairs = int(PHASE3_ASSEMBLY_ROWS[wildcards.sample]["retained_pairs"])
        size = "large" if pairs >= int(PHASE3_ASSEMBLY["large_input_pairs"]) else "small"
        return PHASE3_ASSEMBLY["resources"][f"assembly_{size}"][key]
    return resource



PHASE3_CATALOG_CONFIG_PATH = ROOT / "config" / "phase3_catalog.json"
with PHASE3_CATALOG_CONFIG_PATH.open() as handle:
    PHASE3_CATALOG = json.load(handle)


def catalog_resources(name):
    return {key: PHASE3_CATALOG["resources"][name][key] for key in ("mem_mb", "runtime", "partition")}


PHASE4_CATALOG_CONFIG_PATH = ROOT / "config" / "phase4_catalog.json"
with PHASE4_CATALOG_CONFIG_PATH.open() as handle:
    PHASE4_CATALOG = json.load(handle)


def phase4_resources(name):
    return {key: PHASE4_CATALOG["resources"][name][key] for key in ("mem_mb", "runtime", "partition")}


PHASE5_CONFIG_PATH = ROOT / "config" / "phase5_mapping.json"
with PHASE5_CONFIG_PATH.open() as handle:
    PHASE5 = json.load(handle)


def phase5_resources(name):
    return {key: PHASE5["resources"][name][key] for key in ("mem_mb", "runtime", "partition")}


def phase5_reads(wildcards, mate):
    """Host-depleted cohort input where a host reference exists; regenerated trimmed reads otherwise."""
    if SAMPLES[wildcards.sample]["host_route"] == "none":
        return f"{WORK}/phase5_mapping/trimmed/{wildcards.sample}_R{mate}.fastq.gz"
    return f"{WORK}/phase3_cohort/inputs/{wildcards.sample}_R{mate}.fastq.gz"


def phase5_mapping_resource(key):
    def resource(wildcards):
        size = "large" if SAMPLES[wildcards.sample]["host_route"] == "none" else "small"
        return PHASE5["resources"][f"mapping_{size}"][key]
    return resource


def phase5_scope_samples(wildcards):
    if wildcards.scope == "pilot":
        return PILOT_IDS
    if wildcards.scope == "cohort":
        return SAMPLE_IDS
    raise WorkflowError(f"unknown Phase-5 scope: {wildcards.scope}")


PHASE6_CONFIG_PATH = ROOT / "config" / "phase6_analysis.json"
with PHASE6_CONFIG_PATH.open() as handle:
    PHASE6 = json.load(handle)
with (ROOT / "config" / "phase6_subsets.tsv").open(newline="") as handle:
    PHASE6_SUBSETS = list(csv.DictReader(handle, delimiter="\t"))
PHASE6_SENSITIVITY = sorted({(variant, row["sample_id"]) for row in PHASE6_SUBSETS
                             for variant in row["variants"].split(";")})


def phase6_sensitivity_reads(wildcards, mate):
    """Capped full trimmed reads (shared with Phase 5 trimming) or uncapped trimmed reads."""
    if wildcards.variant == "full_reads":
        return f"{WORK}/phase5_mapping/trimmed/{wildcards.sample}_R{mate}.fastq.gz"
    return f"{WORK}/phase6_analysis/trimmed_uncapped/{wildcards.sample}_R{mate}.fastq.gz"


def phase6_index(wildcards):
    if wildcards.kind == "bacterial":
        base = f"{WORK}/phase4_catalog/{PHASE5['catalog_version']}/bacterial_catalog.fna"
    else:
        base = f"{WORK}/phase5_mapping/{PHASE5['catalog_version']}/viral_catalog.fna"
    return [base + suffix for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")]


EUKARYOTE_CONFIG_PATH = ROOT / "config" / "eukaryote_gate.json"
with EUKARYOTE_CONFIG_PATH.open() as handle:
    EUKARYOTE = json.load(handle)


PILOT_IDS = list(config.get("phase2_pilot_samples", []))
if len(PILOT_IDS) != len(set(PILOT_IDS)):
    raise WorkflowError("phase2 pilot sample IDs are duplicated")
if any(sample not in SAMPLES for sample in PILOT_IDS):
    missing = sorted(set(PILOT_IDS) - set(SAMPLES))
    raise WorkflowError(f"phase2 pilot samples absent from sample table: {missing}")
if config.get("mode") == "production" and not 6 <= len(PILOT_IDS) <= 10:
    raise WorkflowError("production Phase-2 pilot must contain 6--10 samples")


def phase2_validation_paths(sample_ids):
    return expand(f"{WORK}/validation/screens/{{sample}}.json", sample=sample_ids)


def phase2_scope_samples(wildcards):
    if wildcards.scope == "pilot":
        return PILOT_IDS
    if wildcards.scope == "cohort":
        return SAMPLE_IDS
    raise WorkflowError(f"unknown Phase-2 scope: {wildcards.scope}")


_target = config.get("default_target", "phase1_ready")
_target_paths = {
    "phase1_ready": f"{WORK}/stages/phase1.ready",
    "phase1_smoke": f"{WORK}/stages/phase1_smoke.done",
    "screen_pilot": f"{WORK}/stages/screen_pilot.done",
    "screen_cohort": f"{WORK}/stages/screen_cohort.done",
    "phase3_benchmark": f"{WORK}/stages/phase3_benchmark.done",
    "phase3_inputs": f"{WORK}/stages/phase3_inputs.done",
    "phase3_assembly": f"{WORK}/stages/phase3_assembly.done",
    "phase3_catalog": f"{WORK}/stages/phase3_catalog.done",
    "phase4_catalog": f"{WORK}/stages/phase4_catalog.done",
    "phase5_pilot": f"{WORK}/stages/phase5_mapping_pilot.done",
    "phase5_mapping": f"{WORK}/stages/phase5_mapping_cohort.done",
    "phase6_analysis": f"{WORK}/stages/phase6_analysis.done",
    "eukaryote_evidence": f"{WORK}/eukaryote_gate/read_lineages.tsv",
}
if _target not in _target_paths:
    raise WorkflowError(f"unknown default_target: {_target}")
DEFAULT_TARGETS = [_target_paths[_target]]
