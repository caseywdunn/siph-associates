#!/usr/bin/env python3
"""Write stable SHA-256 checksums for Phase-0 generated and source metadata."""
from __future__ import annotations

import hashlib
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUTPUT = REPO / "data" / "metadata" / "manifest.freeze.sha256"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


paths = [
    REPO / "manifest.csv",
    REPO / "data" / "metadata" / "raw_files.tsv",
    REPO / "data" / "metadata" / "library_provenance.tsv",
    REPO / "data" / "metadata" / "input_resources.tsv",
    REPO / "data" / "metadata" / "reference_audits.tsv",
    REPO / "data" / "metadata" / "read_count_corrections.tsv",
    REPO / "data" / "metadata" / "storage_estimate.tsv",
    REPO / "scripts" / "build_manifest.py",
    REPO / "scripts" / "validate_manifest.py",
    REPO / "scripts" / "inventory_phase0_resources.py",
    REPO / "scripts" / "summarize_phase0.py",
    REPO / "scripts" / "freeze_manifest.py",
] + sorted(path for path in (REPO / "data" / "sources").iterdir() if path.is_file())

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text("".join(
    f"{sha256(path)}  {path.relative_to(REPO)}\n" for path in paths
))
print(f"froze {len(paths)} files -> {OUTPUT}")
