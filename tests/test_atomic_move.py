#!/usr/bin/env python3
"""Exercise atomic publication across the node-local and project filesystems."""
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
from common import atomic_move

destination = ROOT / ".cache" / "atomic_move_test.txt"
destination.parent.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(prefix="atomic-move.", dir="/tmp") as temporary:
    source = Path(temporary) / "source.txt"
    source.write_text("phase2-cross-filesystem-test\n")
    atomic_move(source, destination)
    if source.exists() or destination.read_text() != "phase2-cross-filesystem-test\n":
        raise AssertionError("cross-filesystem atomic move failed")
destination.unlink()
print("status\tPASS\n")
