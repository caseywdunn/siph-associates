#!/usr/bin/env python3
"""Verify assembled non-host eukaryotic SSUs against NCBI core_nt with BLAST (URL API).

Each sequence's top hits are recorded with the query date and the database
version reported by NCBI, because NCBI results change as the database grows.
"""
import argparse
import csv
import datetime
import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from xml.etree import ElementTree

URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
parser = argparse.ArgumentParser()
parser.add_argument("--fasta", required=True)
parser.add_argument("--email", required=True)
parser.add_argument("--output", required=True, type=Path)
parser.add_argument("--batch", type=int, default=10)
args = parser.parse_args()


def request(params, data=False):
    encoded = urllib.parse.urlencode({**params, "EMAIL": args.email, "TOOL": "siph_associates"}).encode()
    for attempt in range(6):
        try:
            if data:
                with urllib.request.urlopen(URL, data=encoded, timeout=120) as r:
                    return r.read().decode()
            with urllib.request.urlopen(f"{URL}?{encoded.decode()}", timeout=120) as r:
                return r.read().decode()
        except Exception:
            time.sleep(30 * (attempt + 1))
    raise RuntimeError("NCBI BLAST request failed")


records, name = {}, None
for line in open(args.fasta):
    if line.startswith(">"):
        name = line[1:].split()[0]
        records[name] = []
    elif name:
        records[name].append(line.strip())
names = list(records)
# BLAST parses '|' in query names as database identifiers, so send simple IDs and map back.
alias = {f"q{i:05d}": n for i, n in enumerate(names, start=1)}
local = {n: a for a, n in alias.items()}
rows = []
for start in range(0, len(names), args.batch):
    batch = names[start:start + args.batch]
    query = "".join(f">{local[n]}\n{''.join(records[n])}\n" for n in batch)
    reply = request({"CMD": "Put", "PROGRAM": "blastn", "DATABASE": "core_nt", "QUERY": query,
                     "HITLIST_SIZE": 5}, data=True)
    rid = re.search(r"RID = (\S+)", reply).group(1)
    time.sleep(30)
    while "Status=WAITING" in request({"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}):
        time.sleep(60)
    # BLAST XML is the machine-readable format this API returns reliably (Tabular comes back empty).
    root = ElementTree.fromstring(request({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"}))
    database = root.findtext("BlastOutput_db") or "core_nt"
    seen = set()
    for iteration in root.iter("Iteration"):
        query_id = alias.get(iteration.findtext("Iteration_query-def", "").split()[0], "")
        for rank, hit in enumerate(iteration.iter("Hit"), start=1):
            if rank > 3:
                break
            hsp = hit.find("Hit_hsps/Hsp")
            length = int(hsp.findtext("Hsp_align-len"))
            rows.append({"sequence_id": query_id, "rank": rank, "subject": hit.findtext("Hit_accession"),
                         "title": hit.findtext("Hit_def"),
                         "identity": round(100 * int(hsp.findtext("Hsp_identity")) / length, 2),
                         "alignment_length": length, "evalue": hsp.findtext("Hsp_evalue"),
                         "bitscore": hsp.findtext("Hsp_bit-score"), "database": database,
                         "query_date": datetime.date.today().isoformat(), "rid": rid})
            seen.add(query_id)
    for n in batch:
        if n not in seen:
            rows.append({"sequence_id": n, "rank": 0, "subject": "no significant hit", "title": "", "identity": "",
                         "alignment_length": "", "evalue": "", "bitscore": "", "database": database,
                         "query_date": datetime.date.today().isoformat(), "rid": rid})
    print(f"batch {start // args.batch + 1}: {len(batch)} sequences, RID {rid}", flush=True)

# Titles and organisms for the hit accessions (NCBI nucleotide summaries).
accessions = sorted({r["subject"] for r in rows if r["rank"]})
summaries = {}
for start in range(0, len(accessions), 100):
    reply = urllib.request.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "nuccore", "id": ",".join(accessions[start:start + 100]), "retmode": "json",
         "tool": "siph_associates", "email": args.email}), timeout=120)
    result = json.loads(reply.read().decode())["result"]
    for uid in result.get("uids", []):
        entry = result[uid]
        for key in (entry.get("accessionversion", ""), entry.get("caption", "")):
            summaries[key] = (entry.get("title", ""), entry.get("organism", ""))
    time.sleep(0.4)
for r in rows:
    r["organism"] = summaries.get(r["subject"], ("", ""))[1]
    r["title"] = r["title"] or summaries.get(r["subject"], ("", ""))[0]

fields = ["sequence_id", "rank", "subject", "organism", "title", "identity", "alignment_length", "evalue",
          "bitscore", "database", "query_date", "rid"]
temporary = args.output.with_name(args.output.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
os.replace(temporary, args.output)
