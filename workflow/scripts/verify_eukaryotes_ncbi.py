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
                     "HITLIST_SIZE": 5, "FORMAT_TYPE": "Tabular"}, data=True)
    rid = re.search(r"RID = (\S+)", reply).group(1)
    time.sleep(30)
    while "Status=WAITING" in request({"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}):
        time.sleep(60)
    text = request({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "Tabular", "ALIGNMENTS": 5, "DESCRIPTIONS": 5,
                    "ALIGNMENT_VIEW": "Tabular", "FORMAT_OBJECT": "Alignment"})
    database = re.search(r"# Database: (.+)", text)
    seen = {}
    for line in text.splitlines():
        if line.startswith("#") or "\t" not in line:
            continue
        f = line.split("\t")
        query_id = alias.get(f[0], f[0])
        rank = seen[query_id] = seen.get(query_id, 0) + 1
        if rank <= 3:
            rows.append({"sequence_id": query_id, "rank": rank, "subject": f[1], "identity": f[2],
                         "alignment_length": f[3], "evalue": f[10], "bitscore": f[11].strip(),
                         "database": database.group(1).strip() if database else "core_nt",
                         "query_date": datetime.date.today().isoformat(), "rid": rid})
    for n in batch:
        if n not in seen:
            rows.append({"sequence_id": n, "rank": 0, "subject": "no significant hit", "identity": "",
                         "alignment_length": "", "evalue": "", "bitscore": "",
                         "database": database.group(1).strip() if database else "core_nt",
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
        summaries[entry.get("accessionversion", "")] = (entry.get("title", ""), entry.get("organism", ""))
    time.sleep(0.4)
for r in rows:
    r["title"], r["organism"] = summaries.get(r["subject"], ("", ""))

fields = ["sequence_id", "rank", "subject", "organism", "title", "identity", "alignment_length", "evalue",
          "bitscore", "database", "query_date", "rid"]
temporary = args.output.with_name(args.output.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
os.replace(temporary, args.output)
