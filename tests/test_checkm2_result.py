#!/usr/bin/env python3
"""Check that only CheckM2's exact no-annotation exit becomes an empty result."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
from checkm2_result import NO_ANNOTATIONS_MESSAGE, classify_checkm2_exit

if classify_checkm2_exit(0, "") != "assessed":
    raise AssertionError("successful CheckM2 exit was not classified as assessed")
if classify_checkm2_exit(1, NO_ANNOTATIONS_MESSAGE) != "no_annotations":
    raise AssertionError("exact no-annotation exit was not recognized")
if classify_checkm2_exit(0, NO_ANNOTATIONS_MESSAGE) != "assessed":
    raise AssertionError("a successful exit must not be overridden by stale log text")
try:
    classify_checkm2_exit(1, "ERROR: DIAMOND process crashed")
except RuntimeError:
    pass
else:
    raise AssertionError("an unrelated CheckM2 failure was masked")

print("status\tPASS")
