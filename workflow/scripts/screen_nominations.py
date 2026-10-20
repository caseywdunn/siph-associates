"""Parse completed taxonomic screen tables into the shared nomination schema."""

import csv
import io
import tarfile
from pathlib import Path


def screen_nominations(
    sample_ids: list[str],
    studies: dict[str, str],
    genus_paths: list[str | Path],
    species_paths: list[str | Path],
    sylph_paths: list[str | Path],
    phylo_paths: list[str | Path],
) -> tuple[list[str], list[dict[str, str]]]:
    """Preserve source evidence and reject an incomplete set of per-library inputs."""
    nomination_fields = [
        "sample_id",
        "study",
        "source",
        "rank",
        "candidate_id",
        "candidate_name",
        "estimated_reads",
        "fraction_total_reads",
        "taxonomic_abundance",
        "sequence_abundance",
        "adjusted_ani",
        "true_coverage",
        "raw_support",
    ]
    nominations = []

    for sample, genus_path, species_path, sylph_path, phylo_path in zip(
        sample_ids, genus_paths, species_paths, sylph_paths, phylo_paths, strict=True
    ):
        study = studies[sample]
        for rank, path in (("genus", genus_path), ("species", species_path)):
            with open(path, newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    nominations.append(
                        {
                            "sample_id": sample,
                            "study": study,
                            "source": "Bracken",
                            "rank": rank,
                            "candidate_id": row["taxonomy_id"],
                            "candidate_name": row["name"],
                            "estimated_reads": row["new_est_reads"],
                            "fraction_total_reads": row["fraction_total_reads"],
                            "raw_support": row["kraken_assigned_reads"],
                        }
                    )
        with open(sylph_path, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                nominations.append(
                    {
                        "sample_id": sample,
                        "study": study,
                        "source": "sylph",
                        "rank": "genome",
                        "candidate_id": row["Genome_file"],
                        "candidate_name": Path(row["Genome_file"]).name,
                        "taxonomic_abundance": row["Taxonomic_abundance"],
                        "sequence_abundance": row["Sequence_abundance"],
                        "adjusted_ani": row["Adjusted_ANI"],
                        "true_coverage": row["True_cov"],
                        "raw_support": row["kmers_reassigned"],
                    }
                )
        with tarfile.open(phylo_path, "r:gz") as outer:
            nested_member = next(
                member
                for member in outer.getmembers()
                if member.name.endswith(".phyloFlash.tar.gz")
            )
            nested_stream = outer.extractfile(nested_member)
            if nested_stream is None:
                raise ValueError(f"cannot read nested phyloFlash archive for {sample}")
            with tarfile.open(
                fileobj=io.BytesIO(nested_stream.read()), mode="r:gz"
            ) as inner:
                ntu_member = next(
                    member
                    for member in inner.getmembers()
                    if member.name.endswith(".phyloFlash.NTUfull_abundance.csv")
                )
                ntu_stream = inner.extractfile(ntu_member)
                if ntu_stream is None:
                    raise ValueError(f"cannot read phyloFlash NTU table for {sample}")
                for line in io.TextIOWrapper(ntu_stream, encoding="utf-8"):
                    lineage, count = line.rstrip("\n").rsplit(",", 1)
                    nominations.append(
                        {
                            "sample_id": sample,
                            "study": study,
                            "source": "phyloFlash",
                            "rank": "SSU_lineage",
                            "candidate_id": lineage,
                            "candidate_name": lineage.split(";")[-1],
                            "raw_support": count,
                        }
                    )
    return nomination_fields, nominations
