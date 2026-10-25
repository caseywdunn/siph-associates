#!/usr/bin/env python3
"""Reuse reviewed microbial catalog identities and apply viral interpretation changes."""
from __future__ import annotations
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path
from prepare_coauthor_reporting import apply_classes, read, write

source = Path("figures/catalog_qc/catalog_v3")
out = Path("figures/manuscript_revision/catalog")
out.mkdir(parents=True, exist_ok=True)
inputs = [Path("config/viral_catalog_classifications.tsv"), Path(__file__)]
for name in ["bacterial_catalog.manifest.tsv", "mag_catalog.tsv", "catalog_id_crosswalk.tsv"]:
    src = Path("figures/coauthor_review/catalog") / name
    inputs.append(src)
    shutil.copyfile(src, out / name)
changes = read(inputs[0])
overrides = {r["target_id"]: r for r in changes}
original = []
for name in ["associate", "endogenous_candidate"]:
    src = source / f"viral_{name}.manifest.tsv"
    inputs.append(src)
    original.extend(read(src))
updated = apply_classes(original, overrides, "votu_id")
for name in ["associate", "endogenous_candidate"]:
    write(out / f"viral_{name}.manifest.tsv", [r for r in updated if r["catalog_class"] == name])
summary = {"catalog_counts": dict(Counter(r["catalog_class"] for r in updated)),
           "reclassified_targets": sorted(overrides),
           "sequence_ids_and_checksums_unchanged": True,
           "classification_interpretation": "Possible endogenous candidates; integration is not demonstrated."}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
outputs = sorted(out.glob("*.tsv")) + [out / "summary.json"]
record = {"inputs_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
          "outputs_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs}}
(out / "provenance.json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(summary, indent=2))
