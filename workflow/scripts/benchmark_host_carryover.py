import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, ensure_parents, executable, sha256, version

ensure_parents(list(snakemake.output) + [snakemake.log[0]])
route = str(snakemake.params.route)
if route == "none":
    Path(snakemake.log[0]).write_text("host_route=none; host carryover unavailable without a conspecific reference\n")
    atomic_json(snakemake.output.metrics, {
        "benchmark_id": str(snakemake.wildcards.benchmark_id),
        "sample_id": str(snakemake.params.sample),
        "host_route": route,
        "primary_records": None,
        "host_mapped_primary_records": None,
        "host_mapped_record_fraction": None,
        "reason": "no_conspecific_host_reference",
        "upstream_provenance_sha256": sha256(snakemake.input.upstream),
    })
else:
    prefix = Path(str(snakemake.params.prefix))
    bwa, samtools = executable(prefix, "bwa"), executable(prefix, "samtools")
    environment = os.environ.copy()
    environment["PATH"] = str(prefix / "bin") + os.pathsep + environment.get("PATH", "")
    scratch_parent = Path(snakemake.config["scratch_root"]) / "phase3_benchmark" / "host_qc_tmp"
    scratch_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"hostqc.{snakemake.wildcards.benchmark_id}.",
                                     dir=str(scratch_parent)) as temporary:
        bam = Path(temporary) / "host.bam"
        mapping = [bwa, "mem", "-t", str(snakemake.threads), str(snakemake.params.index_prefix),
                   str(snakemake.input.r1), str(snakemake.input.r2)]
        conversion = [samtools, "view", "-@", str(max(1, snakemake.threads // 2)), "-b", "-o", str(bam), "-"]
        with open(snakemake.log[0], "w") as log:
            mapper = subprocess.Popen(mapping, stdout=subprocess.PIPE, stderr=log, env=environment)
            assert mapper.stdout is not None
            converter = subprocess.Popen(conversion, stdin=mapper.stdout, stdout=log, stderr=log, env=environment)
            mapper.stdout.close()
            converter_code = converter.wait()
            mapper_code = mapper.wait()
        if mapper_code or converter_code:
            raise RuntimeError(f"host carryover mapping failed: bwa={mapper_code}, samtools={converter_code}")
        def count(flags):
            result = subprocess.run([samtools, "view", "-c", "-F", str(flags), str(bam)],
                                    text=True, stdout=subprocess.PIPE, check=True, env=environment)
            return int(result.stdout.strip())
        primary = count(2304)
        mapped = count(2308)
    atomic_json(snakemake.output.metrics, {
        "benchmark_id": str(snakemake.wildcards.benchmark_id),
        "sample_id": str(snakemake.params.sample),
        "host_route": route,
        "primary_records": primary,
        "host_mapped_primary_records": mapped,
        "host_mapped_record_fraction": mapped / primary if primary else 0.0,
        "software": {"bwa": version([bwa]), "samtools": version([samtools, "--version"])},
        "upstream_provenance_sha256": sha256(snakemake.input.upstream),
    })
