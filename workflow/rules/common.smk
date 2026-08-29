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


_target = config.get("default_target", "phase1_ready")
_target_paths = {
    "phase1_ready": f"{WORK}/stages/phase1.ready",
    "phase1_smoke": f"{WORK}/stages/phase1_smoke.done",
    "screen_cohort": f"{WORK}/stages/screen_cohort.done",
}
if _target not in _target_paths:
    raise WorkflowError(f"unknown default_target: {_target}")
DEFAULT_TARGETS = [_target_paths[_target]]
