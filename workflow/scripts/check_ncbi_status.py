#!/usr/bin/env python3
"""Record NCBI Assembly status for GTDB accessions (current and unsuppressed)."""
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


with open(args.genomes, newline="") as handle:
    accessions = sorted(r["genome_id"] for r in csv.DictReader(handle, delimiter="\t") if r["source"] == "GTDB")
records = {}
for start in range(0, len(accessions), 50):
    batch = accessions[start:start + 50]
    term = " OR ".join(f"{a}[Assembly Accession]" for a in batch)
    ids = call("esearch.fcgi", db="assembly", term=term, retmax=500)["esearchresult"]["idlist"]
    if ids:
        result = call("esummary.fcgi", db="assembly", id=",".join(ids))["result"]
        for uid in result["uids"]:
            summary = result[uid]
            records[summary["assemblyaccession"]] = summary

temporary = args.output.with_name(args.output.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["genome_id", "ncbi_found", "ncbi_current", "ncbi_suppressed", "ncbi_species", "ncbi_properties"])
    for accession in accessions:
        summary = records.get(accession)
        if summary is None:
            writer.writerow([accession, "false", "false", "", "", ""])
            continue
        properties = summary.get("propertylist", [])
        suppressed = any("suppressed" in p for p in properties)
        writer.writerow([accession, "true", str("latest" in properties).lower(), str(suppressed).lower(),
                         summary.get("speciesname", ""), ";".join(properties)])
os.replace(temporary, args.output)
print(f"queried {len(accessions)} accessions; found {len(records)}")
