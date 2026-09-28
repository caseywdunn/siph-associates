#!/usr/bin/env python3
"""Screen references, pool with MAG species, draw decoys, and write the frozen catalog."""
import argparse
import csv
import hashlib
import json
import os
from collections import defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
parser.add_argument("--references", required=True, help="nomination table")
parser.add_argument("--decoys", required=True, help="decoy candidate table")
parser.add_argument("--staged", required=True, help="staged genome table")
parser.add_argument("--genomes", required=True, type=Path)
parser.add_argument("--ncbi", required=True)
parser.add_argument("--substitutes", required=True)
parser.add_argument("--substitute-genomes", required=True, type=Path)
parser.add_argument("--checkm2", required=True)
parser.add_argument("--oceandna-gtdbtk", required=True, nargs="+")
parser.add_argument("--mag-catalog", required=True)
parser.add_argument("--mag-representatives", required=True, type=Path)
parser.add_argument("--ani", required=True)
parser.add_argument("--manifest", required=True, type=Path)
parser.add_argument("--fasta", required=True, type=Path)
args = parser.parse_args()

settings = json.loads(Path(args.config).read_text())
qc, derep, decoy_rules = settings["reference_qc"], settings["dereplication"], settings["decoys"]


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


nominations = {r["genome_id"]: r for r in rows(args.references)}
decoys = rows(args.decoys)
staged = {r["genome_id"]: r for r in rows(args.staged)}
ncbi = {r["genome_id"]: r for r in rows(args.ncbi)}
substitutes = {r["genome_id"]: r for r in rows(args.substitutes)}
substitute_min_ani = settings["references"]["substitute_min_ani_to_gtdb_genome"]
quality = {r["Name"]: r for r in rows(args.checkm2)}
placements = {r["user_genome"]: r["classification"] for path in args.oceandna_gtdbtk for r in rows(path)}
domains = {f"d__{d}" for d in qc["oceandna_require_gtdbtk_domains"]}


# Pairwise ANI: species links need both alignment fractions; decoy exclusion uses ANI alone.
species_links, near = defaultdict(set), defaultdict(set)
stem = lambda path: Path(path).name.removesuffix(".fa")
pair_ani = {}
for row in rows(args.ani):
    a, b, ani = stem(row["Ref_file"]), stem(row["Query_file"]), float(row["ANI"])
    pair_ani[frozenset((a, b))] = ani
    if ani >= decoy_rules["exclude_within_ani_of_catalog"]:
        near[a].add(b)
        near[b].add(a)
    if ani >= derep["ani_min"] and min(float(row["Align_fraction_ref"]),
                                       float(row["Align_fraction_query"])) >= derep["align_fraction_min"]:
        species_links[a].add(b)
        species_links[b].add(a)


def screen(genome_id, source):
    """Return (reason or '', completeness, contamination, taxonomy, sequence_id) for a reference or decoy.

    A non-current GTDB record is replaced by its verified current NCBI substitute: the
    substitute must be at least the locked ANI to the GTDB genome and passes the same QC.
    """
    sequence_id = genome_id
    if source == "GTDB":
        status = ncbi.get(genome_id, {})
        taxonomy = nominations.get(genome_id, {}).get("gtdb_taxonomy") or next(
            (d["gtdb_taxonomy"] for d in decoys if d["genome_id"] == genome_id), "")
        if status.get("ncbi_current") != "true":
            substitute = substitutes.get(genome_id, {}).get("substitute_accession")
            if not substitute:
                return "ncbi_not_current_no_substitute", 0.0, 100.0, taxonomy, sequence_id
            if pair_ani.get(frozenset((genome_id, substitute)), 0.0) < substitute_min_ani:
                return "ncbi_substitute_unverified", 0.0, 100.0, taxonomy, sequence_id
            sequence_id = substitute
    else:
        taxonomy = placements.get(genome_id, "not_placed")
        if taxonomy.split(";")[0] not in domains:
            return "no_prokaryotic_placement", 0.0, 100.0, taxonomy, sequence_id
    record = quality.get(sequence_id)
    completeness = float(record["Completeness"]) if record else 0.0
    contamination = float(record["Contamination"]) if record else 100.0
    if record is None:
        return "checkm2_missing", completeness, contamination, taxonomy, sequence_id
    if completeness < qc["min_completeness"] or contamination >= qc["max_contamination_exclusive"]:
        return "checkm2_below_threshold", completeness, contamination, taxonomy, sequence_id
    return "", completeness, contamination, taxonomy, sequence_id


entries = {}
for genome_id, nomination in nominations.items():
    reason, completeness, contamination, taxonomy, sequence_id = screen(genome_id, nomination["source"])
    substitute = substitutes.get(genome_id, {}) if sequence_id != genome_id else {}
    entries[genome_id] = {
        "genome_id": genome_id, "role": "reference", "source": nomination["source"],
        "taxonomy": taxonomy, "completeness": completeness, "contamination": contamination,
        "score": round(completeness - 5 * contamination, 2), "status": reason or "candidate",
        "evidence": f"sylph_libraries={nomination['libraries']};max_ani={nomination['max_ani']};"
                    f"bracken_libraries={nomination['bracken_corroborated_libraries']};"
                    f"phyloflash_libraries={nomination['phyloflash_corroborated_libraries']}",
        "sequence_id": sequence_id,
        "fasta": (args.substitute_genomes if substitute else args.genomes) / f"{sequence_id}.fa",
        "origin": substitute.get("url") or staged[genome_id]["origin"],
        "origin_sha256": substitute.get("fasta_sha256") or staged[genome_id]["origin_sha256"],
        "substitute": f"{sequence_id}:{substitute['substitute_reason']}" if substitute else "",
    }
for mag in rows(args.mag_catalog):
    if mag["catalog_status"] == "retained" and mag["mag_id"] == mag["representative"]:
        entries[mag["species_cluster"]] = {
            "genome_id": mag["species_cluster"], "role": "mag", "source": f"phase3:{mag['mag_id']}",
            "taxonomy": mag["gtdb_classification"], "completeness": float(mag["completeness"]),
            "contamination": float(mag["contamination"]), "score": float(mag["quality_score"]),
            "status": "candidate", "evidence": f"source_library={mag['sample_id']}",
            "sequence_id": mag["species_cluster"], "substitute": "",
            "fasta": args.mag_representatives / f"{mag['species_cluster']}.fa",
            "origin": mag["mag_id"], "origin_sha256": "",
        }

# Greedy centroid dereplication in descending quality order.
pool = sorted((e for e in entries.values() if e["status"] == "candidate"),
              key=lambda e: (-e["score"], e["role"] != "mag", e["genome_id"]))
representatives = []
for entry in pool:
    match = next((r for r in representatives if r["sequence_id"] in species_links[entry["sequence_id"]]), None)
    if match is None:
        representatives.append(entry)
        entry["status"] = "representative"
        entry["members"] = [entry["genome_id"]]
    else:
        entry["status"] = f"clustered_into:{match['genome_id']}"
        match["members"].append(entry["genome_id"])
catalog_genomes = {r["sequence_id"] for r in representatives}

# Decoys: take candidates in seeded order per family until each quota is met.
chosen = []
for family, quota in decoy_rules["families"].items():
    taken = 0
    for candidate in sorted((d for d in decoys if d["family"] == family), key=lambda d: int(d["draw_order"])):
        reason, completeness, contamination, taxonomy, sequence_id = screen(candidate["genome_id"], "GTDB")
        substitute = substitutes.get(candidate["genome_id"], {}) if sequence_id != candidate["genome_id"] else {}
        if not reason and near[sequence_id] & catalog_genomes:
            reason = "within_ani_of_catalog"
        entry = {"genome_id": candidate["genome_id"], "role": "decoy", "source": "GTDB", "taxonomy": taxonomy,
                 "completeness": completeness, "contamination": contamination,
                 "score": round(completeness - 5 * contamination, 2),
                 "evidence": f"decoy_family={family};draw_order={candidate['draw_order']}",
                 "sequence_id": sequence_id,
                 "fasta": (args.substitute_genomes if substitute else args.genomes) / f"{sequence_id}.fa",
                 "origin": substitute.get("url") or staged[candidate["genome_id"]]["origin"],
                 "origin_sha256": substitute.get("fasta_sha256") or staged[candidate["genome_id"]]["origin_sha256"],
                 "substitute": f"{sequence_id}:{substitute['substitute_reason']}" if substitute else ""}
        if taken >= quota:
            continue
        if reason:
            entry["status"] = f"decoy_excluded:{reason}"
        else:
            entry["status"] = "representative"
            entry["members"] = [entry["genome_id"]]
            taken += 1
        entries[f"decoy:{candidate['genome_id']}"] = entry
        chosen.append(entry)
    if taken < quota:
        raise SystemExit(f"decoy family {family} filled {taken} of {quota}")

# Stable catalog IDs: real genomes by taxonomy, then decoys; contigs prefixed by catalog ID.
members = sorted(representatives, key=lambda e: (e["taxonomy"], e["genome_id"]))
decoy_members = sorted((e for e in chosen if e["status"] == "representative"),
                       key=lambda e: (e["taxonomy"], e["genome_id"]))
for i, entry in enumerate(members, start=1):
    entry["catalog_id"] = f"BAC{i:05d}"
for i, entry in enumerate(decoy_members, start=1):
    entry["catalog_id"] = f"DEC{i:03d}"

temporary = args.fasta.with_name(args.fasta.name + ".tmp")
with temporary.open("w") as out:
    for entry in members + decoy_members:
        bases = contigs = 0
        for line in open(entry["fasta"]):
            if line.startswith(">"):
                contigs += 1
                out.write(f">{entry['catalog_id']}|{line[1:].split()[0]}\n")
            else:
                bases += len(line.strip())
                out.write(line)
        entry["contigs"], entry["bases"] = contigs, bases
        entry["fasta_sha256"] = sha256(entry["fasta"])
os.replace(temporary, args.fasta)

fields = ["catalog_id", "genome_id", "role", "source", "status", "substitute", "taxonomy", "completeness",
          "contamination", "score", "contigs", "bases", "cluster_members", "evidence", "origin", "origin_sha256",
          "fasta_sha256"]
temporary = args.manifest.with_name(args.manifest.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    ordered = members + decoy_members + sorted(
        (e for e in entries.values() if e["status"] != "representative"), key=lambda e: (e["role"], e["genome_id"]))
    for entry in ordered:
        writer.writerow({**entry, "catalog_id": entry.get("catalog_id", ""),
                         "cluster_members": ",".join(entry.get("members", [])),
                         "contigs": entry.get("contigs", ""), "bases": entry.get("bases", ""),
                         "fasta_sha256": entry.get("fasta_sha256", "")})
os.replace(temporary, args.manifest)
print(f"catalog genomes={len(members)} (mags={sum(e['role'] == 'mag' for e in members)}, "
      f"references={sum(e['role'] == 'reference' for e in members)}) decoys={len(decoy_members)}")
