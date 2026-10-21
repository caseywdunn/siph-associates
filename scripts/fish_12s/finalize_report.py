#!/usr/bin/env python3
"""Validate and combine the completed fish screens without searching reads again."""

import hashlib
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from functools import cache
from pathlib import Path

from Bio import SeqIO

from pilot import digest as file_digest
from summarize_cohort import table, write_table

digest = cache(file_digest)


def verify_report(path: Path) -> dict:
    record = json.loads((path / "run.json").read_text())
    for filename, expected in record["output_sha256"].items():
        assert digest(path / filename) == expected, (path, filename)
    for filename, expected in record["input_sha256"].items():
        assert digest(Path(filename)) == expected, filename
    return record


def main() -> None:
    out = Path("figures/fish_12s")
    out.mkdir(parents=True, exist_ok=True)
    review = Path("data/results/post_expansion/2026-10-21/reporting/fish_12s")
    review.mkdir(parents=True, exist_ok=True)
    host = Path("data/results/fish_12s/final_review/2026-10-19/host_depleted")
    full = Path("data/results/fish_12s/full_reads/2026-10-18/review")
    libraries = table(Path("config/fish_12s_full_reads.tsv"))
    directories = [host] + [full / r["analysis_sample_id"] for r in libraries]
    assert len(libraries) == 44
    additions = table(Path("data/metadata/manko2026/sample_metadata.tsv"))
    added_root = Path("data/results/bargmannia/2026-10-19/fish_12s/review")
    directories += [added_root / r["analysis_sample_id"] for r in additions]
    assert len(additions) == 8
    samples, pairs, sources = [], [], {}
    for directory in directories:
        verify_report(directory)
        sources[str(directory / "run.json")] = digest(directory / "run.json")
        mode = "host_depleted" if directory == host else "capped_full_reads"
        samples.extend(
            dict(r, read_source=mode) for r in table(directory / "sample_summary.tsv")
        )
        pairs.extend(table(directory / "pair_review.tsv"))
    expected = {
        r["analysis_sample_id"]
        for r in table(
            Path(
                "data/results/post_expansion/2026-10-21/reporting/metadata/sample_metadata.tsv"
            )
        )
    }
    assert len(samples) == len({r["analysis_sample_id"] for r in samples}) == 213
    assert {r["analysis_sample_id"] for r in samples} == expected
    assert len(pairs) == len(
        {(r["analysis_sample_id"], r["read_pair_id"]) for r in pairs}
    )
    counts = Counter((r["analysis_sample_id"], r["review_category"]) for r in pairs)
    for sample in samples:
        sid = sample["analysis_sample_id"]
        categories = [
            "long_close_fish_match",
            "other_fish_similarity",
            "other_or_no_hit",
        ]
        assert sum(int(sample[k]) for k in categories) == int(sample["candidate_pairs"])
        for category in categories:
            assert int(sample[category]) == counts[sid, category], (sid, category)

    positives = [
        dict(r) for r in pairs if r["review_category"] == "long_close_fish_match"
    ]
    wanted = {(r["analysis_sample_id"], r["read_pair_id"]) for r in positives}
    sequences = defaultdict(dict)
    config = json.loads(Path("config/fish_12s_host_depleted_cohort.json").read_text())
    reuse = {r["analysis_sample_id"] for r in table(Path(config["reuse_manifest"]))}
    for sid in sorted({s for s, _ in wanted}):
        relative = Path(sid, "host_depleted")
        if sid.startswith("Manko2026__"):
            folders = [added_root.parent / "union" / sid / "capped_trimmed"]
        elif sid in {r["analysis_sample_id"] for r in libraries}:
            folders = [full.parent / "union" / sid / "capped_trimmed"]
        elif sid in reuse:
            folders = [Path(config["reuse_screen"]) / relative]
            extra = (
                Path("data/results/fish_12s/host_depleted_aligner_only/2026-10-17")
                / relative
            )
            if extra.exists():
                folders.append(extra)
        else:
            folders = [Path(config["run"]) / "union" / relative]
        for folder in folders:
            mapping = folder / (
                "candidate_reads.tsv" if sid in reuse else "candidate_pairs.tsv"
            )
            queries = {}
            mate_counts = Counter()
            for row in table(mapping):
                rid = row["read_pair_id"]
                if sid in reuse:
                    mate_counts[rid] += 1
                    queries[row["query_id"]] = (
                        rid,
                        int(row.get("mate", mate_counts[rid])),
                    )
                else:
                    for mate in (1, 2):
                        queries[f"{row['candidate_id']}_{mate}"] = (rid, mate)
            fasta = folder / "candidates.fasta"
            sources[str(fasta)] = digest(fasta)
            for record in SeqIO.parse(fasta, "fasta"):
                rid, mate = queries[record.id]
                if (sid, rid) in wanted:
                    sequences[sid, rid][mate] = str(record.seq).upper()
    seq_groups = defaultdict(list)
    for row in positives:
        key = row["analysis_sample_id"], row["read_pair_id"]
        seq = sequences[key]
        assert set(seq) == {1, 2}, key
        # Treat swapped mates as the same paired sequence; no UMI inference.
        signature = hashlib.sha256("|".join(sorted(seq.values())).encode()).hexdigest()
        seq_groups[signature].append(key)
        row["paired_sequence_sha256"] = signature
    for row in positives:
        group = seq_groups[row["paired_sequence_sha256"]]
        row["identical_pairs_in_cohort"] = len(group)
        row["flowcell"] = row["read_pair_id"].split(":")[2]

    summary = {
        "libraries": len(samples),
        "host_depleted_libraries": sum(
            r["read_source"] == "host_depleted" for r in samples
        ),
        "full_read_libraries": len(libraries) + len(additions),
        "pairs_screened": sum(int(r["pairs_screened"]) for r in samples),
        "full_read_pairs_screened": sum(
            int(r["pairs_screened"])
            for r in samples
            if r["read_source"] == "capped_full_reads"
        ),
        "categories": dict(Counter(r["review_category"] for r in pairs)),
        "positive_libraries": len({r["analysis_sample_id"] for r in positives}),
        "unique_positive_paired_sequences": len(seq_groups),
        "identical_positive_pair_groups": sum(len(v) > 1 for v in seq_groups.values()),
        "method_recovery": {
            m: sum(int(r[m]) for r in positives)
            for m in ("regex", "bbduk25", "bowtie2", "sortmerna")
        },
        "interpretation": "Long, close competitive mitochondrial matches; read-pair counts are not prey individuals or independent molecules. Exact paired-sequence duplication checked; PCR duplicates with sequence errors and sequencing contamination are not excluded. Species identities remain reference-match hypotheses.",
    }
    write_table(out / "sample_summary.tsv", samples)
    write_table(out / "positive_pairs.tsv", positives)
    write_table(review / "all_candidate_pairs.tsv", pairs)
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    values = {
        "Libraries": len(samples),
        "PositiveLibraries": summary["positive_libraries"],
        "StrongPairs": len(positives),
        "FullReadLibraries": len(libraries) + len(additions),
        "FullReadPairs": summary["full_read_pairs_screened"],
    }
    (out / "values.tex").write_text(
        "".join(f"\\newcommand{{\\fish{k}}}{{{v:,}}}\n" for k, v in values.items())
    )
    sources[__file__] = digest(Path(__file__))
    (out / "provenance.json").write_text(
        json.dumps(
            {
                "completed_utc": datetime.now(timezone.utc).isoformat(),
                "command": sys.argv,
                "inputs_sha256": sources,
                "outputs_sha256": {
                    p.name: digest(p)
                    for p in out.iterdir()
                    if p.name != "provenance.json"
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
