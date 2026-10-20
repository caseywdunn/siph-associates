#!/usr/bin/env python3
"""Extend frozen nominations with completed Bargmannia screens without rereading FASTQs."""

import json
import sys
from pathlib import Path

from assembly_completion import sha256, table, write_table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
from screen_nominations import screen_nominations


def main() -> None:
    root = Path("data/results/bargmannia/2026-10-19")
    output = Path("data/results/post_expansion/2026-10-20/aggregation/cohort")
    samples = [r["analysis_sample_id"] for r in table("config/bargmannia_samples.tsv")]
    inputs = []
    for sample in samples:
        trim = root / "provenance/trim" / f"{sample}.json"
        for method in ("kraken_bracken", "sylph", "phyloflash"):
            path = root / "provenance" / method / f"{sample}.json"
            record = json.loads(path.read_text())
            if record["sample_id"] != sample or record[
                "upstream_provenance_sha256"
            ] != sha256(trim):
                raise ValueError(
                    f"Screen provenance does not match trim input: {sample} {method}"
                )
            products = record.get(
                "outputs", [record.get("output", record.get("archive"))]
            )
            for item in products:
                if (
                    method == "kraken_bracken"
                    and item
                    and "/bracken/" not in item["path"]
                ):
                    continue
                if not item or Path(item["path"]).stat().st_size != item["size_bytes"]:
                    raise ValueError(
                        f"Screen product differs from completed producer: {sample} {method}"
                    )
                inputs.append(Path(item["path"]))
            inputs.append(path)

    def paths(kind, suffix):
        return [str(root / "screens" / kind / f"{s}{suffix}") for s in samples]

    fields, new = screen_nominations(
        samples,
        dict.fromkeys(samples, "Manko2026"),
        paths("bracken", ".G.tsv"),
        paths("bracken", ".S.tsv"),
        paths("sylph", ".profile.tsv"),
        paths("phyloflash", ".tar.gz"),
    )
    source = Path("data/results/aggregation/cohort/candidate_nominations.tsv")
    expected = json.loads(Path("config/phase4_catalog.json").read_text())["nomination"][
        "source_sha256"
    ]
    if sha256(source) != expected:
        raise ValueError("Original nomination table differs from its frozen identity")
    original = table(source)
    if any(r["sample_id"] in samples for r in original):
        raise ValueError("New library already included in original nominations")
    write_table(output / "candidate_nominations.tsv", original + new, fields)
    (output / "provenance.json").write_text(
        json.dumps(
            {
                "original_libraries": 205,
                "new_libraries": samples,
                "new_nomination_rows": len(new),
                "inputs_sha256": {
                    str(p): sha256(p) for p in sorted(set(inputs + [source]))
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(f"Added {len(new)} nomination rows from {len(samples)} libraries")


if __name__ == "__main__":
    main()
