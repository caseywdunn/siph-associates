#!/usr/bin/env python3
"""Export reviewed supporting data without rerunning analyses.

The selection pins original analysis paths, sizes and SHA-256 checksums.
Use --dry-run before exporting, or --check to verify a public clone without
requiring the original analysis outputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELECTION = Path("config/supporting_tables.tsv")
OUT = Path("tables")
MAX_FILE_BYTES = 5 * 1024 * 1024
MAX_TOTAL_BYTES = 20 * 1024 * 1024


def sha256(path: Path) -> str:
    """Return a file checksum without loading the whole file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def retained_path(relative: str) -> Path:
    """Reject selections that leave the analysis repository."""
    path = (ROOT / relative).resolve()
    if Path(relative).is_absolute() or not path.is_relative_to(ROOT):
        raise ValueError(f"Path leaves the repository: {relative}")
    return path


def read_selection() -> list[dict]:
    """Read the bounded, explicitly reviewed file selection."""
    with (ROOT / SELECTION).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    destinations = set()
    for row in rows:
        row["bytes"] = int(row["bytes"])
        retained_path(row["source"])
        target = retained_path(row["destination"])
        if not target.is_relative_to(ROOT / OUT):
            raise ValueError(f"Destination is outside tables/: {target}")
        if row["destination"] in destinations:
            raise ValueError(f"Duplicate destination: {target}")
        destinations.add(row["destination"])
        if not 0 < row["bytes"] <= MAX_FILE_BYTES:
            raise ValueError(f"File exceeds the selection size limit: {target}")
    if not rows or sum(row["bytes"] for row in rows) > MAX_TOTAL_BYTES:
        raise ValueError("Selection is empty or exceeds the total size limit")
    return rows


def verify(rows: list[dict], field: str) -> None:
    """Check every file before writing any exported data."""
    for row in rows:
        path = retained_path(row[field])
        if path.stat().st_size != row["bytes"] or sha256(path) != row["sha256"]:
            raise ValueError(f"Reviewed size or checksum does not match: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true", help="verify original files")
    mode.add_argument("--check", action="store_true", help="verify exported files")
    args = parser.parse_args()
    rows = read_selection()
    verify(rows, "destination" if args.check else "source")
    if args.check:
        provenance = json.loads((ROOT / OUT / "provenance.json").read_text())
        if provenance["selection_sha256"] != sha256(ROOT / SELECTION):
            raise ValueError("Export provenance does not match the selection")
        if provenance["artifacts"] != rows:
            raise ValueError("Export provenance does not match the retained files")
    elif not args.dry_run:
        revision = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
        dirty = bool(
            subprocess.check_output(
                ["git", "status", "--porcelain", "--untracked-files=no"],
                cwd=ROOT,
                text=True,
            ).strip()
        )
        provenance = {
            "schema_version": 1,
            "source_base_git_revision": revision,
            "source_worktree_dirty": dirty,
            "exported_utc": datetime.now(timezone.utc).isoformat(),
            "command": "python3 scripts/export_supporting_tables.py",
            "selection": str(SELECTION),
            "selection_sha256": sha256(ROOT / SELECTION),
            "script_sha256": sha256(Path(__file__)),
            "artifacts": rows,
        }
        with tempfile.TemporaryDirectory(prefix="supporting-tables-", dir=ROOT) as tmp:
            stage = Path(tmp)
            for row in rows:
                target = stage / row["destination"]
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(retained_path(row["source"]), target)
                if sha256(target) != row["sha256"]:
                    raise ValueError(f"Staged copy checksum failed: {target}")
            for row in rows:
                target = retained_path(row["destination"])
                target.parent.mkdir(parents=True, exist_ok=True)
                os.replace(stage / row["destination"], target)
            manifest = stage / "provenance.json"
            manifest.write_text(json.dumps(provenance, indent=2) + "\n")
            os.replace(manifest, ROOT / OUT / "provenance.json")
        verify(rows, "destination")
    print(
        json.dumps(
            {
                "mode": "check"
                if args.check
                else "dry-run"
                if args.dry_run
                else "export",
                "files": len(rows),
                "bytes": sum(row["bytes"] for row in rows),
                "status": "verified",
            }
        )
    )


if __name__ == "__main__":
    main()
