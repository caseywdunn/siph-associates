#!/usr/bin/env python3
"""Require passed COI verification and the exact retained input inventory."""

import argparse
import csv
import hashlib
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    gate = Path(config["verification_run"]) / "gate.json"
    checks = json.loads(gate.read_text())
    if not checks["passed"] or not all(checks["checks"].values()):
        raise ValueError("COI scientific verification has not passed")
    with Path(config["manifest"]).open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 161 or len({r["analysis_sample_id"] for r in rows}) != 161:
        raise ValueError("Unexpected cohort")
    for row in rows:
        metrics = json.loads(Path(row["upstream_metrics"]).read_text())
        if metrics["strategy"] != "reference_depleted" or metrics[
            "retained_pairs"
        ] != int(row["pairs"]):
            raise ValueError(f"Changed input counts for {row['analysis_sample_id']}")
        for mate in [1, 2]:
            path = Path(row[f"r{mate}_path"])
            if path.stat().st_size != int(row[f"r{mate}_bytes"]):
                raise ValueError(f"Changed input file: {path}")
    out = Path(config["run"])
    out.mkdir(parents=True, exist_ok=True)
    (out / "VERIFICATION_PASSED").write_text(
        json.dumps(
            dict(
                gate_sha256=hashlib.sha256(gate.read_bytes()).hexdigest(),
                manifest_sha256=hashlib.sha256(
                    Path(config["manifest"]).read_bytes()
                ).hexdigest(),
                libraries=len(rows),
                pairs=sum(int(r["pairs"]) for r in rows),
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
