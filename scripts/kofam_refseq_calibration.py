#!/usr/bin/env python3
"""Calibrate the relaxed KofamScan tier against curated RefSeq annotation.

Four complete RefSeq Metamycoplasmataceae genomes from the gene-content set are
compared gene by gene (pyrodigal CDS matched to RefSeq CDS by contig, strand and
stop coordinate). For each relaxed score fraction, the script reports recall of
curated focal genes (by RefSeq product name) and, as a precision proxy, how often
an assigned KO's definition shares a word with the RefSeq product.

The single Mollicutes 2-oxoacid dehydrogenase E1 scores best against the
branched-chain profiles, so K00166/K00167 count as matches for PDH E1 alpha/beta.

Reads data/results/mycoplasmatales_gene_content/ and downloads RefSeq GFFs from
NCBI Datasets. Writes data/results/mycoplasmatales_gene_content/calibration/.

Usage: python3 scripts/kofam_refseq_calibration.py
"""
from __future__ import annotations

import csv
import io
import re
import urllib.request
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENE = ROOT / "data" / "results" / "mycoplasmatales_gene_content"
OUT = GENE / "calibration"
GENOMES = {
    "GCF_000183385.1": "Mycoplasmopsis bovis",
    "GCF_000085865.1": "Metamycoplasma hominis",
    "GCF_900660575.1": "Mycoplasmopsis pulmonis",
    "GCF_000008365.1": "Mycoplasma mobile",
}
EVALUE_MAX = 1e-5
FRACTIONS = (None, 0.75, 0.6, 0.5)  # None = strict thresholds only
# Curated focal genes: RefSeq product pattern -> accepted KOs.
CURATED = {
    "PDH E1 alpha": (r"pyruvate dehydrogenase.*alpha|alpha-ketoacid dehydrogenase subunit alpha", {"K00161", "K00166"}),
    "PDH E1 beta": (r"pyruvate dehydrogenase.*beta|alpha-ketoacid dehydrogenase subunit beta", {"K00162", "K00167"}),
    "PDH E2": (r"dihydrolipoamide acetyltransferase|dihydrolipoyllysine.*acetyltransferase", {"K00627"}),
    "PDH E3": (r"dihydrolipoyl dehydrogenase", {"K00382"}),
    "pta": (r"phosphate acetyltransferase", {"K00625"}),
    "ackA": (r"^acetate kinase", {"K00925"}),
    "arcA": (r"arginine deiminase", {"K01478"}),
    "arcB": (r"ornithine carbamoyltransferase", {"K00611"}),
    "arcC": (r"carbamate kinase", {"K00926"}),
    "glpK": (r"^glycerol kinase", {"K00864"}),
}
STOP = {"protein", "subunit", "family", "domain", "containing", "putative", "type", "system", "component", "chain",
        "enzyme", "the", "and", "like", "factor", "dependent", "related", "specific"}


def words(text):
    return {w for w in re.findall(r"[a-z0-9]{4,}", text.lower()) if w not in STOP}


def refseq_products(accession):
    """RefSeq CDS products keyed by (contig, strand, stop coordinate), cached under OUT/refseq/."""
    gff = OUT / "refseq" / f"{accession}.gff"
    if not gff.exists():
        gff.parent.mkdir(parents=True, exist_ok=True)
        url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"
               "?include_annotation_type=GENOME_GFF")
        with zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(url).read())) as archive:
            name = next(n for n in archive.namelist() if n.endswith(".gff"))
            gff.write_bytes(archive.read(name))
    products = {}
    for line in open(gff):
        fields = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(fields) < 9 or fields[2] != "CDS":
            continue
        match = re.search(r"product=([^;]+)", fields[8])
        stop = fields[4] if fields[6] == "+" else fields[3]
        products[(fields[0], fields[6], stop)] = urllib.request.unquote(match.group(1)) if match else ""
    return products


hits = defaultdict(list)
for line in open(GENE / "kofamscan.tsv"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if f[1].split("|")[0].removeprefix("GTDB_") in GENOMES:
        hits[f[1]].append((f[0] == "*", f[2], float(f[4]), float(f[3]) if f[3] else None, float(f[5])))
definitions = {}
for line in open(GENE / "kofam" / "ko_list"):
    f = line.rstrip("\n").split("\t")
    definitions[f[0]] = f[-1]


def assign(protein, fraction):
    candidates = hits.get(protein, [])
    passing = [h for h in candidates if h[0]]
    if passing:
        return max(passing, key=lambda h: h[2])[1], "strict"
    if fraction is None:
        return None, None
    near = [h for h in candidates if h[3] and h[4] <= EVALUE_MAX and h[2] >= fraction * h[3]]
    return (max(near, key=lambda h: h[2])[1], "relaxed") if near else (None, None)


genes = []  # (genome, protein id, RefSeq product or None)
for accession in GENOMES:
    products = refseq_products(accession)
    index = 0
    for line in open(GENE / "genes" / f"GTDB_{accession}.gff"):
        fields = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(fields) < 9 or fields[2] != "CDS":
            continue
        index += 1
        stop = fields[4] if fields[6] == "+" else fields[3]
        genes.append((accession, f"GTDB_{accession}|{fields[0]}_{index}", products.get((fields[0], fields[6], stop))))

summary, per_gene = [], []
for fraction in FRACTIONS:
    label = "strict" if fraction is None else f"relaxed_{fraction:g}"
    found, total, agree, assigned = Counter(), Counter(), Counter(), Counter()
    for accession, protein, product in genes:
        if product is None:
            continue
        ko, tier = assign(protein, fraction)
        for gene, (pattern, accepted) in CURATED.items():
            if re.search(pattern, product, re.I):
                total[gene] += 1
                found[gene] += ko in accepted
                per_gene.append({"rule": label, "genome": accession, "protein": protein, "curated_gene": gene,
                                 "refseq_product": product, "assigned_ko": ko or "", "tier": tier or "",
                                 "recovered": ko in accepted})
        if ko and "hypothetical" not in product.lower():
            assigned[tier] += 1
            agree[tier] += bool(words(definitions.get(ko, "")) & words(product))
    summary.append({
        "rule": label, "curated_recovered": sum(found.values()), "curated_total": sum(total.values()),
        "missed": ";".join(f"{g}:{total[g] - found[g]}" for g in CURATED if total[g] - found[g]),
        "strict_assignments": assigned["strict"],
        "strict_name_agreement": round(agree["strict"] / assigned["strict"], 3) if assigned["strict"] else "",
        "relaxed_added_assignments": assigned["relaxed"],
        "relaxed_name_agreement": round(agree["relaxed"] / assigned["relaxed"], 3) if assigned["relaxed"] else "",
    })

OUT.mkdir(parents=True, exist_ok=True)
for name, records in (("recall_by_rule.tsv", summary), ("curated_genes.tsv", per_gene)):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
matched = sum(p is not None for _, _, p in genes)
print(f"pyrodigal genes {len(genes)}, matched to RefSeq CDS {matched}")
for row in summary:
    print("\t".join(str(row[k]) for k in row))
