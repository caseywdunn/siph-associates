#!/usr/bin/env python3
"""Exercise immutable input staging between project storage and scratch."""
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
from common import stage_readonly_input

source = ROOT / ".cache" / "stage_readonly_input.source.txt"
source.parent.mkdir(parents=True, exist_ok=True)
source.write_text("phase3-cross-filesystem-test\n")
try:
    with tempfile.TemporaryDirectory(prefix="stage-input.", dir="/tmp") as temporary:
        destination = Path(temporary) / "input.txt"
        stage_readonly_input(source, destination)
        if destination.read_text() != source.read_text():
            raise AssertionError("staged input does not match its source")
        if source.stat().st_dev != destination.parent.stat().st_dev and not destination.is_symlink():
            raise AssertionError("cross-filesystem input was not staged as a symbolic link")
finally:
    source.unlink(missing_ok=True)

print("status\tPASS")
