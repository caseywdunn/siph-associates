#!/usr/bin/env python3
"""Summarize competitive cross-domain matches without forcing a domain assignment."""

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from prepare import domain, read_parents, write_tsv

HIT_FIELDS = [
    "query_id",
    "subject_id",
    "pident",
    "length",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
    "title",
]


def rows(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def classify_hits(
    hits: list[dict], fraction: float, minimum_aa: int, maximum_evalue: float
) -> tuple[str, list[dict]]:
    usable = [
        h
        for h in hits
        if float(h["evalue"]) <= maximum_evalue and int(h["length"]) >= minimum_aa
    ]
    if not usable:
        return "no_qualifying_hit", []
    best = max(float(h["bitscore"]) for h in usable)
    near = [h for h in usable if float(h["bitscore"]) >= fraction * best]
    labels = set().union(*(h["domains"] for h in near))
    if "unresolved_taxonomy" in labels:
        return "unresolved_taxonomy", near
    if len(labels) != 1:
        return "ambiguous", near
    return next(iter(labels)), near


def annotate_subjects(hits: list[dict], reference: Path, source: str) -> None:
    if source == "swissprot":
        wanted = {h["subject_id"] for h in hits}
        metadata = {}
        with (reference / "swissprot_metadata.tsv").open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                for key in [row["subject_id"], row["accession"]]:
                    if key in wanted:
                        metadata[key] = row
        for hit in hits:
            record = metadata[hit["subject_id"]]
            hit["domains"] = {record["domain"]}
            hit["organisms"] = record["organism"]
        return
    # nr was built without accession-to-taxid data. Resolve bracketed organism
    # names against NCBI taxonomy, retaining missing/ambiguous names explicitly.
    wanted = set()
    for hit in hits:
        names = re.findall(r"\[([^\[\]]+)\]", hit["title"])
        hit["organisms"] = ";".join(names)
        wanted.update(names)
    parents = read_parents(reference / "taxonomy/nodes.dmp")
    name_taxids = defaultdict(set)
    with (reference / "taxonomy/names.dmp").open() as handle:
        for line in handle:
            fields = [f.strip() for f in line.split("|")]
            if fields[1] in wanted:
                name_taxids[fields[1]].add(int(fields[0]))
    for hit in hits:
        labels = set()
        for name in hit["organisms"].split(";"):
            taxids = name_taxids.get(name)
            labels.update(domain(t, parents) for t in taxids) if taxids else labels.add(
                "unresolved_taxonomy"
            )
        hit["domains"] = labels


def report(
    config: dict,
    query_path: Path,
    hit_path: Path,
    source: str,
    out: Path,
    controls: bool,
) -> None:
    queries = rows(query_path)
    query_ids = {r["query_id"] for r in queries}
    if len(query_ids) != len(queries):
        raise ValueError("Duplicate query metadata")
    with hit_path.open() as handle:
        hits = list(csv.DictReader(handle, fieldnames=HIT_FIELDS, delimiter="\t"))
    if any(h["query_id"] not in query_ids or None in h for h in hits):
        raise ValueError("Alignment columns or query IDs disagree with metadata")
    root = Path(config["run"])
    annotate_subjects(
        hits, Path(config.get("reference_root", root / "reference")), source
    )
    by_query = defaultdict(list)
    for hit in hits:
        by_query[hit["query_id"]].append(hit)
    results = []
    for query in queries:
        label, near = classify_hits(
            by_query[query["query_id"]],
            config["near_best_fraction"],
            config["minimum_alignment_aa"],
            config["maximum_evalue"],
        )
        best = max(near, key=lambda h: float(h["bitscore"])) if near else None
        results.append(
            {
                **query,
                "assignment": label,
                "near_best_hits": len(near),
                "near_best_domains": ";".join(
                    sorted(set().union(*(h["domains"] for h in near)))
                ),
                "near_best_organisms": ";".join(sorted({h["organisms"] for h in near})),
                "best_subject": best["subject_id"] if best else "",
                "best_title": best["title"] if best else "",
                "best_pident": best["pident"] if best else "",
                "best_alignment_aa": best["length"] if best else "",
                "best_bitscore": best["bitscore"] if best else "",
                "best_evalue": best["evalue"] if best else "",
            }
        )
    write_tsv(out, results, list(results[0]))
    counts = Counter((r["bin_id"], r["assignment"]) for r in results)
    write_tsv(
        out.with_suffix(".counts.tsv"),
        [
            {"bin_id": b, "assignment": d, "queries": n}
            for (b, d), n in sorted(counts.items())
        ],
        ["bin_id", "assignment", "queries"],
    )
    if controls:
        clear = Counter()
        for query in results:
            expected = query["expected_domain"]
            near = by_query[query["query_id"]]
            exact = [
                h
                for h in near
                if float(h["pident"]) == 100
                and int(h["length"]) == int(query["length"])
            ]
            if not exact or not any(expected in h["domains"] for h in exact):
                raise ValueError(
                    f"Control lacks its full-length exact reference: {query['query_id']}"
                )
            if query["assignment"] not in {expected, "ambiguous"}:
                raise ValueError(f"Control domain is incorrect: {query}")
            clear[expected] += query["assignment"] == expected
        if any(clear[d] < 1 for d in ["Bacteria", "Archaea", "Eukaryota"]):
            raise ValueError(f"No unambiguous control for a domain: {clear}")
        (out.parent / "CONTROLS_PASSED").write_text(
            json.dumps(dict(clear), indent=2) + "\n"
        )
    elif query_path.name == "genes.tsv":
        assignments = {r["query_id"]: r for r in results}
        markers = [
            {**m, **assignments[m["query_id"]]}
            for m in rows(root / "archaeal_marker_hits.tsv")
        ]
        if markers:
            write_tsv(
                out.with_suffix(".archaeal_markers.tsv"), markers, list(markers[0])
            )
        contig_counts = Counter(
            (r["bin_id"], r["contig_id"], r["assignment"]) for r in results
        )
        write_tsv(
            out.with_suffix(".contigs.tsv"),
            [
                {"bin_id": b, "contig_id": c, "assignment": d, "genes": n}
                for (b, c, d), n in sorted(contig_counts.items())
            ],
            ["bin_id", "contig_id", "assignment", "genes"],
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--queries", required=True, type=Path)
    parser.add_argument("--hits", required=True, type=Path)
    parser.add_argument("--source", choices=["swissprot", "nr"], required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--controls", action="store_true")
    args = parser.parse_args()
    report(
        json.loads(args.config.read_text()),
        args.queries,
        args.hits,
        args.source,
        args.output,
        args.controls,
    )


if __name__ == "__main__":
    main()
