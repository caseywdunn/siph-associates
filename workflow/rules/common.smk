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
}
if _target not in _target_paths:
    raise WorkflowError(f"unknown default_target: {_target}")
DEFAULT_TARGETS = [_target_paths[_target]]
