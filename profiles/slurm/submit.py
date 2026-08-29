#!/usr/bin/env python3
"""Submit one Snakemake jobscript with declared rule resources."""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

from snakemake.utils import read_job_properties

jobscript = sys.argv[-1]
properties = read_job_properties(jobscript)
resources = properties.get("resources", {})
threads = int(properties.get("threads", 1))
rule = re.sub(r"[^A-Za-z0-9_.-]", "_", properties.get("rule", "job"))[:60]
wildcards = properties.get("wildcards", {})
sample = re.sub(r"[^A-Za-z0-9_.-]", "_", str(wildcards.get("sample", wildcards.get("route", ""))))[:60]
name = f"sa.{rule}" + (f".{sample}" if sample else "")
log_dir = Path("logs") / "slurm"
log_dir.mkdir(parents=True, exist_ok=True)

command = [
    "sbatch", "--parsable", f"--job-name={name}",
    f"--partition={resources.get('partition', 'day')}",
    f"--cpus-per-task={threads}",
    f"--mem={int(resources.get('mem_mb', 4000))}M",
    f"--time={int(resources.get('runtime', 60))}",
    f"--output={log_dir}/{name}.%j.out",
    jobscript,
]
result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, check=False)
if result.returncode:
    sys.stderr.write(result.stderr)
    raise SystemExit(result.returncode)
print(result.stdout.strip().split(";")[0])
