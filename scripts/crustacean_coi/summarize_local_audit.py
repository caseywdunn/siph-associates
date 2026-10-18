#!/usr/bin/env python3
"""Document the two NA35 false-positive COI calls using bacterial catalog matches."""

import argparse
import csv
import re
from pathlib import Path


def table(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict]) -> None:
    with path.open("w") as handle:
        writer = csv.DictWriter(
            handle, list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    args = parser.parse_args()
    metadata = {r["catalog_id"]: r for r in table(args.catalog)}
    selected = {"c00000008": "Alteromonas", "c00000016": "Vibrio"}
    matches = []
    with (args.root / "audit/bacterial_catalog.sam").open() as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            row = line.rstrip().split("\t")
            flag = int(row[1])
            if flag & (4 | 256 | 2048) or not any(
                f"NA35__{cid}_" in row[0] for cid in selected
            ):
                continue
            tags = {t.split(":", 2)[0]: t.split(":", 2)[2] for t in row[11:]}
            aligned = sum(
                int(n)
                for n, op in re.findall(r"(\d+)([MIDNSHP=X])", row[5])
                if op in "M=XID"
            )
            reference = metadata[row[2].split("|")[0]]
            matches.append(
                dict(
                    query=row[0],
                    reference=row[2],
                    genome_id=reference["genome_id"],
                    taxonomy=reference["taxonomy"],
                    strand="-" if flag & 16 else "+",
                    position=int(row[3]),
                    mapping_quality=int(row[4]),
                    cigar=row[5],
                    aligned_nt=aligned,
                    edit_distance=int(tags["NM"]),
                    identity_percent=round(
                        100 * (aligned - int(tags["NM"])) / aligned, 3
                    ),
                    alignment_score=int(tags["AS"]),
                )
            )
    assert len(matches) == 4, "Expected primary alignments for both mates of two pairs"
    write(args.root / "audit/na35_bacterial_matches.tsv", matches)
    original = {
        r["candidate_id"]: r
        for r in table(args.root / "reclassified/Ahuja2026__NA35/classified_pairs.tsv")
    }
    pairs = []
    for cid, organism in selected.items():
        evidence = [r for r in matches if f"__{cid}_" in r["query"]]
        assert len(evidence) == 2 and len({r["reference"] for r in evidence}) == 1
        left, right = sorted(evidence, key=lambda r: r["position"])
        assert left["strand"] == "+" and right["strand"] == "-"
        assert right["position"] - left["position"] < 1000
        assert all(
            r["identity_percent"] >= 90
            and r["aligned_nt"] >= 130
            and r["mapping_quality"] >= 30
            and f"g__{organism};" in r["taxonomy"]
            for r in evidence
        )
        pairs.append(
            dict(
                sample="Ahuja2026__NA35",
                candidate_id=cid,
                read_pair_id=original[cid]["read_pair_id"],
                nucleotide_panel_label=original[cid]["classification"],
                audited_label="non_target_bacterial",
                evidence=f"Both mates map in compatible orientation on the same {organism} catalog contig; see na35_bacterial_matches.tsv",
                interpretation="False positive for crustacean prey; no species-level bacterial claim",
            )
        )
    write(args.root / "audit/na35_audited_pairs.tsv", pairs)


if __name__ == "__main__":
    main()
