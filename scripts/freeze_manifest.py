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


def freeze_input_paths(repo: Path) -> list[Path]:
    """Include enrichment inputs and implementation, excluding historical freezes."""
    relative_paths = [
        "manifest.csv",
        "data/metadata/raw_files.tsv",
        "data/metadata/library_provenance.tsv",
        "data/metadata/input_resources.tsv",
        "data/metadata/reference_audits.tsv",
        "data/metadata/read_count_corrections.tsv",
        "data/metadata/storage_estimate.tsv",
        "scripts/build_manifest.py",
        "scripts/validate_manifest.py",
        "scripts/inventory_phase0_resources.py",
        "scripts/summarize_phase0.py",
        "scripts/freeze_manifest.py",
        "scripts/merge_sample_metadata.py",
    ]
    paths = {repo / relative for relative in relative_paths}
    paths.update(path for path in (repo / "data/sources").iterdir() if path.is_file())
    enrichment = repo / "data/sources/metadata_enrichment"
    if enrichment.exists():
        paths.update(path for path in enrichment.rglob("*") if path.is_file())
    for relative in (
        "data/metadata/sample_metadata_updates.tsv",
        "scripts/extract_sample_sheet_metadata.py",
        "scripts/review_sample_metadata.py",
        "scripts/score_physalia_depths.py",
        "docs/metadata_depth_conventions.md",
    ):
        path = repo / relative
        if path.exists():
            paths.add(path)
    return sorted(paths)


def main() -> None:
    paths = freeze_input_paths(REPO)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        "".join(f"{sha256(path)}  {path.relative_to(REPO)}\n" for path in paths)
    )
    print(f"froze {len(paths)} files -> {OUTPUT}")


if __name__ == "__main__":
    main()
