#!/usr/bin/env python3
"""Prepare/check a two-lane fish screen through the production cap/QC route."""

import argparse
import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/fish_12s"))
from full_read_inputs import table, write_table
from pilot import records

ROOT = Path("data/results/fish_12s/full_reads_validation/2026-10-18")
FISH = {
    "motif_R1",
    "motif_reverse_R2",
    "fish_away_from_motif",
    "cartilaginous_motif_variant",
}


def prepare() -> None:
    source = table(Path("config/fish_12s_cohort_validation.tsv"))[0]
    originals = [list(records(Path(source[f"r{mate}_path"]))) for mate in [1, 2]]
    ROOT.mkdir(parents=True, exist_ok=True)
    paths = {1: [], 2: []}
    for lane in [1, 2]:
        for mate in [1, 2]:
            path = ROOT / f"lane{lane}_R{mate}.fastq.gz"
            paths[mate].append(str(path.resolve()))
            with gzip.open(path, "wt") as handle:
                for rep in range(5):
                    for title, seq, qual in originals[mate - 1]:
                        name = title.split("/")[0]
                        handle.write(
                            f"@{name}:lane{lane}:rep{rep}/{mate}\n{seq}\n+\n{qual}\n"
                        )
    for label, cap in [("capped", 60), ("uncapped", 200)]:
        manifest = ROOT / f"{label}.tsv"
        write_table(
            manifest,
            [
                {
                    "analysis_sample_id": label,
                    "host": "synthetic",
                    "raw_pairs": 90,
                    "r1_paths": ";".join(paths[1]),
                    "r2_paths": ";".join(paths[2]),
                    "previous_trim_provenance": "",
                }
            ],
        )
        config = json.loads(Path("config/fish_12s_full_reads.json").read_text())
        config.update(libraries=str(manifest), run=str(ROOT / label), cap_pairs=cap)
        (ROOT / f"{label}.json").write_text(json.dumps(config, indent=2) + "\n")


def check() -> None:
    for label, cap in [("capped", 60), ("uncapped", 90)]:
        directory = ROOT / label
        metrics = json.loads((directory / f"trim/{label}.json").read_text())
        assert metrics["pairs_entering_fastp"] == cap
        assert metrics["pairs_after_fastp"] <= cap
        candidates = table(
            directory / f"union/{label}/capped_trimmed/candidate_pairs.tsv"
        )
        reviewed = table(directory / f"review/{label}/pair_review.tsv")
        assert len(reviewed) == len(candidates)
        assert {r["read_pair_id"].split(":")[0] for r in reviewed} == FISH | {"frog"}
        for row in reviewed:
            control = row["read_pair_id"].split(":")[0]
            expected = "long_close_fish_match" if control in FISH else "other_or_no_hit"
            assert row["review_category"] == expected, row
        print(label, cap, metrics["pairs_after_fastp"], len(candidates), "PASS")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["prepare", "check"])
    args = parser.parse_args()
    prepare() if args.action == "prepare" else check()
