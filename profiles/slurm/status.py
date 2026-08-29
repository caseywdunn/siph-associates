#!/usr/bin/env python3
"""Return Snakemake's generic-cluster status for one SLURM job."""
from __future__ import annotations

import subprocess
import sys

job_id = sys.argv[1]
result = subprocess.run(
    ["sacct", "-n", "-P", "-j", job_id, "--format=State"],
    text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False,
)
states = [line.split("|")[0].split()[0].split("+")[0] for line in result.stdout.splitlines() if line.strip()]
if not states:
    queued = subprocess.run(["squeue", "-h", "-j", job_id], stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, check=False)
    print("running" if queued.stdout else "running")
elif any(state in {"FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "BOOT_FAIL"} for state in states):
    print("failed")
elif all(state == "COMPLETED" for state in states):
    print("success")
else:
    print("running")
