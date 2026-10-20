#!/usr/bin/env python3
"""Describe completed fish-screen alignments without rerunning sequence searches.

These are review bins, not validated presence or species assignments. A long,
close fish match spans >=100 nt at >=95% identity, with only fish references
within 5% of the best bit score and no truncated near-best hit list. RefSeq fish
membership uses accession-linked GenBank lineage annotations; all reference
names and competing evidence remain available for manual review.
Whole mitochondrial references can extend beyond the recruited 12S fragment.
"""

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from classify_candidates import evidence, pair_category, read_hits, reference_metadata
from pilot import digest, read_id, records


def table(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    with path.open("w") as handle:
        writer = csv.DictWriter(
            handle, list(rows[0]) if rows else fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("config/fish_12s_host_depleted_cohort.json")
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument(
        "--reference-metadata",
        type=Path,
        default=Path(
            "data/results/fish_12s/reference_metadata/2026-10-18/references.tsv"
        ),
    )
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    manifest_path = args.manifest or Path(config["manifest"])
    manifest = table(manifest_path)
    reuse = (
        {r["analysis_sample_id"] for r in table(Path(config["reuse_manifest"]))}
        if config.get("reuse_manifest")
        else set()
    )
    meta = reference_metadata(args.reference_metadata)
    inputs = {
        args.config,
        manifest_path,
        args.reference_metadata,
        Path(__file__).with_name("classify_candidates.py"),
    }
    if reuse:
        inputs.add(Path(config["reuse_manifest"]))
    pairs, samples = [], []
    for sample in manifest:
        sid = sample["analysis_sample_id"]
        relative = Path(sid, sample.get("read_pair_id", "host_depleted"))
        directory = (
            Path(config["reuse_screen"]) / relative
            if sid in reuse
            else Path(config["run"]) / "union" / relative
        )
        directories = [directory]
        extra = (
            Path("data/results/fish_12s/host_depleted_aligner_only/2026-10-17")
            / relative
        )
        if sid in reuse and extra.exists():
            directories.append(extra)
        paired = defaultdict(dict)
        methods = defaultdict(lambda: dict(regex=0, bbduk25=0, bowtie2=0, sortmerna=0))
        for directory in directories:
            run = directory / "run.json"
            assert json.loads(run.read_text())["status"] == "complete", run
            inputs.add(run)
            mapping = directory / (
                "candidate_reads.tsv" if sid in reuse else "candidate_pairs.tsv"
            )
            inputs.add(mapping)
            queries = {}
            extra_mates = Counter()
            for row in table(mapping):
                rid = row["read_pair_id"]
                if sid in reuse:
                    extra_mates[rid] += 1
                    mate = int(row.get("mate", extra_mates[rid]))
                    queries[row["query_id"]] = (rid, mate)
                    for key, source in [
                        ("regex", "regex"),
                        ("bbduk25", "bbduk25_pair"),
                    ]:
                        methods[rid][key] |= int(row.get(source, 0))
                else:
                    for mate in [1, 2]:
                        queries[f"{row['candidate_id']}_{mate}"] = (rid, mate)
                    methods[rid].update({key: int(row[key]) for key in methods[rid]})
            path = directory / "alignment/hits.tsv"
            inputs.add(path)
            hits = read_hits(path, meta)
            assert not (set(hits) - queries.keys()), (path, "unknown query")
            for qid, (rid, mate) in queries.items():
                assert mate in [1, 2] and mate not in paired[rid], (sid, rid, mate)
                paired[rid][mate] = evidence(hits.get(qid, []))
        if sid in reuse:
            aligners = Path(config["reuse_aligners"]) / relative
            for path in [aligners / "sortmerna.fq", aligners / "bowtie2.sam"]:
                inputs.add(path)
            for title, _, _ in records(aligners / "sortmerna.fq"):
                assert read_id(title) in paired
                methods[read_id(title)]["sortmerna"] = 1
            with (aligners / "bowtie2.sam").open() as handle:
                for line in handle:
                    if not line.startswith("@"):
                        fields = line.split("\t")
                        if not int(fields[1]) & 4:
                            assert read_id(fields[0]) in paired
                            methods[read_id(fields[0])]["bowtie2"] = 1
        counts, close_methods = Counter(), Counter()
        for rid, mates in paired.items():
            assert set(mates) == {1, 2}, (sid, rid)
            category = pair_category(mates)
            if category == "long_close_fish_match":
                close_methods.update(methods[rid])
            counts[category] += 1
            row = dict(
                analysis_sample_id=sid,
                read_pair_id=rid,
                review_category=category,
                **methods[rid],
            )
            for mate, e in sorted(mates.items()):
                row[f"R{mate}_category"] = e["category"]
                for key in [
                    "reference_id",
                    "name",
                    "identity",
                    "length",
                    "query_length",
                    "subject_start",
                    "subject_end",
                    "marker_overlap_nt",
                ]:
                    row[f"R{mate}_{key}"] = e["top"].get(key, "")
                row[f"R{mate}_near_best_names"] = ";".join(e["near_names"])
                row[f"R{mate}_hit_list_truncated"] = e["truncated"]
            pairs.append(row)
        samples.append(
            dict(
                analysis_sample_id=sid,
                host=sample["host"],
                pairs_screened=int(sample["pairs"]),
                candidate_pairs=len(paired),
                **{
                    key: counts[key]
                    for key in [
                        "long_close_fish_match",
                        "other_fish_similarity",
                        "other_or_no_hit",
                    ]
                },
                **{
                    f"long_close_{method}": close_methods[method]
                    for method in ["regex", "bbduk25", "bowtie2", "sortmerna"]
                },
            )
        )
    assert (
        len(samples) == len({r["analysis_sample_id"] for r in samples}) == len(manifest)
    )
    args.output.mkdir(parents=True, exist_ok=True)
    write_table(args.output / "sample_summary.tsv", samples)
    fields = [
        "analysis_sample_id",
        "read_pair_id",
        "review_category",
        "regex",
        "bbduk25",
        "bowtie2",
        "sortmerna",
    ]
    for mate in [1, 2]:
        fields.extend(
            f"R{mate}_{key}"
            for key in [
                "category",
                "reference_id",
                "name",
                "identity",
                "length",
                "query_length",
                "subject_start",
                "subject_end",
                "marker_overlap_nt",
                "near_best_names",
                "hit_list_truncated",
            ]
        )
    write_table(args.output / "pair_review.tsv", pairs, fields)
    provenance = dict(
        completed_utc=datetime.now(timezone.utc).isoformat(),
        command=sys.argv,
        python=sys.version,
        interpretation=__doc__,
        source_sha256=digest(Path(__file__)),
        input_sha256={str(p): digest(p) for p in sorted(inputs)},
        output_sha256={
            name: digest(args.output / name)
            for name in ["sample_summary.tsv", "pair_review.tsv"]
        },
    )
    (args.output / "run.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(
        json.dumps(
            dict(
                libraries=len(samples),
                pairs_screened=sum(s["pairs_screened"] for s in samples),
                categories=Counter(p["review_category"] for p in pairs),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
