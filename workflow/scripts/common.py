"""Shared, standard-library-only helpers for atomic workflow jobs."""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path


def ensure_parents(paths):
    for path in paths:
        Path(str(path)).parent.mkdir(parents=True, exist_ok=True)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, payload):
    path = Path(str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def normalize_read_id(header):
    token = header.strip().split()[0]
    if token.startswith("@"):
        token = token[1:]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def fastq_pair_count(r1, r2):
    counts = []
    for path in (r1, r2):
        lines = 0
        with gzip.open(path, "rt", errors="strict") as handle:
            for lines, _ in enumerate(handle, start=1):
                pass
        if lines % 4:
            raise ValueError(f"FASTQ has incomplete record: {path} ({lines} lines)")
        counts.append(lines // 4)
    if counts[0] != counts[1]:
        raise ValueError(f"mate counts differ: {r1}={counts[0]}, {r2}={counts[1]}")
    return counts[0]


def executable(prefix, name):
    path = Path(str(prefix)) / "bin" / name
    if not path.is_file():
        raise FileNotFoundError(f"missing executable: {path}")
    return str(path)


def resolve_executable(prefix, name):
    """Prefer an activated pinned environment, then the declared cluster prefix."""
    active = shutil.which(name)
    return active if active else executable(prefix, name)


def version(command):
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, check=False)
    return result.stdout.strip().splitlines()[0] if result.stdout.strip() else "unknown"


def file_record(path, checksum=False):
    path = Path(str(path))
    record = {"path": str(path), "size_bytes": path.stat().st_size}
    if checksum:
        record["sha256"] = sha256(path)
    return record
