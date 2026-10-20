#!/usr/bin/env python3
"""Shared competitive-alignment review criteria for controls and real samples."""

import csv
from collections import defaultdict
from pathlib import Path


def reference_metadata(path: Path) -> dict[str, dict]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {r["reference_id"]: r for r in rows}
    if len(result) != len(rows):
        raise ValueError(f"Duplicate reference identifiers: {path}")
    return result


def read_hits(path: Path, references: dict[str, dict]) -> dict[str, list[dict]]:
    hits = defaultdict(list)
    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            qid, rid = row[:2]
            lookup = rid.split("|")[1] if rid.startswith("ref|") else rid
            annotation = references[
                lookup
            ]  # Unknown references must not be silently classified.
            lo, hi = sorted([int(row[8]), int(row[9])])
            intervals = [
                tuple(map(int, s.split("-")))
                for s in annotation["marker_intervals_1based"].split(";")
                if s
            ]
            overlap = sum(
                max(0, min(hi, end) - max(lo, start) + 1) for start, end in intervals
            )
            hits[qid].append(
                dict(
                    reference_id=rid,
                    name=annotation["name"],
                    fish=annotation["fish"] == "1",
                    identity=float(row[2]),
                    length=int(row[3]),
                    bitscore=float(row[11]),
                    query_length=int(row[12]),
                    subject_start=int(row[8]),
                    subject_end=int(row[9]),
                    marker_overlap_nt=overlap,
                )
            )
    return hits


def evidence(hits: list[dict]) -> dict:
    if not hits:
        return dict(category="no_hit", top={}, near_names=[], truncated=False)
    hits = sorted(hits, key=lambda h: h["bitscore"], reverse=True)
    best = hits[0]
    near = [h for h in hits if h["bitscore"] >= best["bitscore"] * 0.95]
    truncated = len(hits) >= 200 and hits[-1]["bitscore"] >= best["bitscore"] * 0.95
    close = best["identity"] >= 95 and best["length"] >= 100
    if all(h["fish"] for h in near) and close and not truncated:
        category = "long_close_fish_match"
    elif any(h["fish"] for h in near):
        category = "other_fish_similarity"
    else:
        category = "other_reference"
    return dict(
        category=category,
        top=best,
        near_names=sorted({h["name"] for h in near}),
        truncated=truncated,
    )


def pair_category(mates: dict[int, dict]) -> str:
    categories = {m["category"] for m in mates.values()}
    conflicting = any(
        m["category"] == "other_reference"
        and m["top"]["identity"] >= 95
        and m["top"]["length"] >= 100
        for m in mates.values()
    )
    if "long_close_fish_match" in categories and not conflicting:
        return "long_close_fish_match"
    if categories & {"long_close_fish_match", "other_fish_similarity"}:
        return "other_fish_similarity"
    return "other_or_no_hit"
