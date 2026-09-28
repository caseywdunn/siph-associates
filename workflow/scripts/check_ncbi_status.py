#!/usr/bin/env python3
"""Record NCBI Assembly status for GTDB accessions and name current substitutes.

An accession is judged on its own database: a GenBank (GCA) record must be
latest_genbank, a RefSeq (GCF) record latest_refseq and not suppressed_refseq.
A non-current record is substituted by the identical current record in the other
database, or by the latest version of the same accession. The substitute's
genome is verified against the GTDB genome downstream.
"""
import argparse
import csv
import json
import os
import time
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
parser = argparse.ArgumentParser()
parser.add_argument("--genomes", required=True, help="staged genome table")
parser.add_argument("--email", required=True)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()


def call(endpoint, **params):
    params.update(retmode="json", tool="siph_associates", email=args.email)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}", timeout=60) as r:
                time.sleep(0.4)  # stay under three requests per second without an API key
                return json.load(r)
        except Exception:
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"NCBI request failed: {endpoint} {params}")


def summaries(terms):
    """Assembly docsums for search terms, keyed by every accession each record carries."""
    found = {}
    for start in range(0, len(terms), 50):
        batch = terms[start:start + 50]
        ids = call("esearch.fcgi", db="assembly", term=" OR ".join(f"{t}[Assembly Accession]" for t in batch),
                   retmax=1000)["esearchresult"]["idlist"]
        for chunk in range(0, len(ids), 200):
            result = call("esummary.fcgi", db="assembly", id=",".join(ids[chunk:chunk + 200]))["result"]
            for uid in result["uids"]:
                record = result[uid]
                synonyms = record.get("synonym", {})
                for accession in {record["assemblyaccession"], synonyms.get("genbank"), synonyms.get("refseq")}:
                    if accession:
                        found[accession] = record
    return found


def side(accession):
    return "refseq" if accession.startswith("GCF") else "genbank"


def current(record, accession):
    properties = record.get("propertylist", [])
    return f"latest_{side(accession)}" in properties and f"suppressed_{side(accession)}" not in properties


def ftp_path(record, accession):
    return record.get(f"ftppath_{side(accession)}", "").replace("ftp://", "https://")


with open(args.genomes, newline="") as handle:
    accessions = sorted(r["genome_id"] for r in csv.DictReader(handle, delimiter="\t") if r["source"] == "GTDB")
records = summaries(accessions)
# Newer versions of non-current accessions: search the unversioned accession.
stale = [a for a in accessions if a not in records or not current(records[a], a)]
# Search both databases' unversioned accessions (GCA_x and GCF_x share digits).
bases = {a.rsplit(".", 1)[0] for a in stale}
bases |= {("GCA" if b.startswith("GCF") else "GCF") + b[3:] for b in bases}
versions = summaries(sorted(bases)) if stale else {}

temporary = args.output.with_name(args.output.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["genome_id", "ncbi_found", "ncbi_current", "ncbi_suppressed", "ncbi_species",
                     "substitute_accession", "substitute_reason", "substitute_ftp", "ncbi_properties"])
    for accession in accessions:
        record = records.get(accession)
        properties = record.get("propertylist", []) if record else []
        is_current = bool(record) and current(record, accession)
        substitute, reason, ftp = "", "", ""
        if record and not is_current:
            other = record.get("synonym", {}).get("refseq" if side(accession) == "genbank" else "genbank")
            if other and record.get("synonym", {}).get("similarity") == "identical" and current(record, other):
                substitute, reason, ftp = other, "identical_other_database", ftp_path(record, other)
        if not is_current and not substitute:
            digits = accession.rsplit(".", 1)[0][3:]
            newer = sorted((a for a in versions if a[3:].rsplit(".", 1)[0] == digits and a != accession
                            and current(versions[a], a)),
                           key=lambda a: (int(a.rsplit(".", 1)[1]), a.startswith("GCF")))
            if newer:
                substitute, reason = newer[-1], "newer_version"
                ftp = ftp_path(versions[substitute], substitute)
        writer.writerow([accession, str(bool(record)).lower(), str(is_current).lower(),
                         str(f"suppressed_{side(accession)}" in properties).lower(),
                         record.get("speciesname", "") if record else "", substitute, reason, ftp,
                         ";".join(properties)])
os.replace(temporary, args.output)
print(f"queried {len(accessions)} accessions; not current {len(stale)}")
