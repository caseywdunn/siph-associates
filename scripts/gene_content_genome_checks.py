#!/usr/bin/env python3
"""Per-contig checks of the siphonophore Mycoplasmatales genomes in the gene-content set.

For each contig: length, GC, coding density (pyrodigal, code 4), genes, and genes with an
in-frame UGA. Mycoplasmatales read UGA as Trp, so their genes carry internal UGA codons,
while genes from a code-11 contaminant essentially never do. The GTDB Mycoplasmopsis bovis
genome is included as a reference.

Reads data/results/mycoplasmatales_gene_content/{genomes,genes}; writes
data/results/mycoplasmatales_gene_content/checks/contig_checks.tsv and genome_checks.tsv.

Usage: python3 scripts/gene_content_genome_checks.py
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENE = ROOT / "data" / "results" / "mycoplasmatales_gene_content"
OUT = GENE / "checks"
GENOMES = ["MAGSP0005", "MAGSP0007", "MAGSP0010", "MAGSP0011", "MAGSP0012", "MAGSP0029", "MAGSP0031",
           "GTDB_GCF_000183385.1"]
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def write(path, records):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


contig_rows, genome_rows = [], []
for genome in GENOMES:
    sequences, name = {}, None
    for line in open(GENE / "genomes" / f"{genome}.fa"):
        if line.startswith(">"):
            name = line[1:].split()[0]
            sequences[name] = []
        else:
            sequences[name].append(line.strip().upper())
    sequences = {k: "".join(v) for k, v in sequences.items()}
    coding, genes, uga = defaultdict(int), defaultdict(int), defaultdict(int)
    for line in open(GENE / "genes" / f"{genome}.gff"):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[2] != "CDS":
            continue
        start, end = int(f[3]), int(f[4])
        gene = sequences[f[0]][start - 1:end]
        if f[6] == "-":
            gene = gene.translate(COMPLEMENT)[::-1]
        codons = [gene[i:i + 3] for i in range(3, len(gene) - 3, 3)]  # internal codons only
        coding[f[0]] += end - start + 1
        genes[f[0]] += 1
        uga[f[0]] += "TGA" in codons
    for contig, sequence in sequences.items():
        contig_rows.append({"genome": genome, "contig": contig, "length": len(sequence),
                            "gc": round((sequence.count("G") + sequence.count("C")) / len(sequence), 4),
                            "coding_density": round(coding[contig] / len(sequence), 4),
                            "genes": genes[contig], "genes_with_inframe_uga": uga[contig]})
    rows = [r for r in contig_rows if r["genome"] == genome]
    big = [r for r in rows if r["genes"] >= 10]
    genome_rows.append({
        "genome": genome, "contigs": len(rows), "genes": sum(genes.values()),
        "fraction_genes_with_inframe_uga": round(sum(uga.values()) / sum(genes.values()), 3),
        "contigs_ge10_genes_without_uga": sum(r["genes_with_inframe_uga"] == 0 for r in big),
        "contig_gc_min": min(r["gc"] for r in rows), "contig_gc_max": max(r["gc"] for r in rows),
        "coding_density": round(sum(coding.values()) / sum(len(s) for s in sequences.values()), 3),
    })

OUT.mkdir(parents=True, exist_ok=True)
write(OUT / "contig_checks.tsv", contig_rows)
write(OUT / "genome_checks.tsv", genome_rows)
for row in genome_rows:
    print("\t".join(str(v) for v in row.values()))
