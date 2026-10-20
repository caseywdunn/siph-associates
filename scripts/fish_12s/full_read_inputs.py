#!/usr/bin/env python3
"""Freeze reference-free libraries, or inventory one regenerated paired input."""

import argparse
import csv
import json
from pathlib import Path


def table(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(
            handle, list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def freeze(output: Path) -> None:
    metadata = {
        r["analysis_sample_id"]: r
        for r in table(Path("figures/cohort/sample_metadata.tsv"))
    }
    completed = {
        r["analysis_sample_id"]
        for r in table(Path("config/prey_marker_host_depleted.tsv"))
    }
    rows = []
    for row in table(Path("config/samples.tsv")):
        if row["host_route"] != "none":
            continue
        sid = row["sample_id"]
        meta = metadata[sid]
        if sid in completed or int(row["read_pairs"]) != int(meta["read_pairs"]):
            raise ValueError(f"Inconsistent cohort membership/count: {sid}")
        previous = Path(f"data/results/provenance/phase5_mapping/trim/{sid}.json")
        provenance = json.loads(previous.read_text())
        ordered_paths = row["r1_paths"].split(";") + row["r2_paths"].split(";")
        if ordered_paths != [p["path"] for p in provenance["inputs"]]:
            raise ValueError(f"Lane order differs from the original seeded cap: {sid}")
        paths = {p["path"]: p for p in provenance["inputs"]}
        for mate in [1, 2]:
            for path in row[f"r{mate}_paths"].split(";"):
                if Path(path).stat().st_size != paths[path]["size_bytes"]:
                    raise ValueError(f"Raw file differs from earlier inventory: {path}")
        rows.append(
            dict(
                analysis_sample_id=sid,
                host=meta["species_current"],
                raw_pairs=row["read_pairs"],
                r1_paths=row["r1_paths"],
                r2_paths=row["r2_paths"],
                previous_trim_provenance=str(previous),
            )
        )
    if (
        len(rows) != 44
        or len({r["analysis_sample_id"] for r in rows} | completed) != 205
    ):
        raise ValueError("Expected 44 additional and 205 total libraries")
    priority = {
        f"Ahuja2024__{sid}": i
        for i, sid in enumerate(["CWD146", "CWD131", "CWD20", "CWD11"])
    }
    rows.sort(
        key=lambda r: (
            priority.get(r["analysis_sample_id"], 4),
            r["analysis_sample_id"],
        )
    )
    write_table(output, rows)
    print(
        f"{len(rows)} libraries; {sum(min(int(r['raw_pairs']), 200_000_000) for r in rows):,} capped raw pairs"
    )


def prepared(libraries: Path, sample: str, provenance_path: Path, output: Path) -> None:
    row = next(r for r in table(libraries) if r["analysis_sample_id"] == sample)
    provenance = json.loads(provenance_path.read_text())
    if provenance["sample_id"] != sample or provenance["raw_pairs_expected"] != int(
        row["raw_pairs"]
    ):
        raise ValueError(
            "Regenerated read provenance does not match the selected library"
        )
    if row.get("previous_trim_provenance"):
        previous = json.loads(Path(row["previous_trim_provenance"]).read_text())
        for key in ["parameters", "pairs_entering_fastp", "pairs_after_fastp"]:
            if previous[key] != provenance[key]:
                raise ValueError(
                    f"Regenerated read preparation differs from original: {sample} {key}"
                )
        if {p["path"]: p["sha256"] for p in previous["inputs"]} != {
            p["path"]: p["sha256"] for p in provenance["inputs"]
        }:
            raise ValueError(f"Raw input identity changed: {sample}")
    result = dict(
        analysis_sample_id=sample,
        host=row["host"],
        pairs=provenance["pairs_after_fastp"],
        read_source="capped_quality_filtered",
        read_pair_id="capped_trimmed",
        host_route="none",
        upstream_metrics=str(provenance_path),
    )
    for mate, record in enumerate(provenance["outputs"], 1):
        path = Path(record["path"])
        if path.stat().st_size != record["size_bytes"]:
            raise ValueError(f"Prepared input changed: {path}")
        result[f"r{mate}_path"] = str(path.resolve())
        result[f"r{mate}_bytes"] = record["size_bytes"]
    write_table(output, [result])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--libraries", type=Path)
    parser.add_argument("--sample")
    parser.add_argument("--provenance", type=Path)
    args = parser.parse_args()
    if args.sample:
        prepared(args.libraries, args.sample, args.provenance, args.output)
    else:
        freeze(args.output)


if __name__ == "__main__":
    main()
