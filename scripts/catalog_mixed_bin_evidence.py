#!/usr/bin/env python3
"""Evidence that v1 catalog genome BAC00033 was a mixed bin, and a screen of v2 MAGs for the same problem.

BAC00033 (v1) merged into the complete clade-A genome BAC00025 in v2. Its contigs are split by
whether they align (minimap2 asm20, >=50% of length) to BAC00025. For each part the script reports
coding density (Prodigal), alignment to the P. physalis host assembly and to the v2 catalog, and
reads per library from the v1 BAMs. "Orphan" libraries are those where v1 called BAC00033 present
but v2 calls BAC00025 none or trace.

For every v2 MAG representative, contigs <50% coding are measured, aligned to the host assembly,
and their share of reads across the v2 BAMs is reported. Presence needs >=10% breadth, so such
contigs alone cannot create a presence call when they are a smaller share of the genome.

Needs about 48 GB of memory (host-genome index). Usage, from the project root:
  sbatch -p day -c 8 --mem=48G -t 4:00:00 --wrap "python3 scripts/catalog_mixed_bin_evidence.py"
Writes data/results/phase4_catalog/v2/bin_evidence/.
"""
from __future__ import annotations

import csv
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "data" / "results"
OUT = RESULTS / "phase4_catalog" / "v2" / "bin_evidence"
TOOLS = Path("/gpfs/gibbs/project/dunn/cwd7/conda_envs/assembly_env/bin")
HOST = Path("/gpfs/ycga/work/dunn/cwd7/databases/refgenomes/physalia_physalis/"
            "GCA_041430235.2_Physalia_TX2017-38_primary_02_genomic.fna")
MIXED, MERGED_INTO = "BAC00033", "BAC00025"
PRESENT = {"validated", "high_confidence"}


def rows(path, delimiter="\t"):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write(name, records):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def read_fasta(path, keep):
    sequences, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            name = name if keep(name) else None
            if name:
                sequences[name] = []
        elif name:
            sequences[name].append(line.strip())
    return {k: "".join(v) for k, v in sequences.items()}


def write_fasta(path, sequences):
    with open(path, "w") as out:
        for name, sequence in sequences.items():
            out.write(f">{name}\n{sequence}\n")


def aligned_fraction(query_fasta, target, tmp, preset="asm20"):
    """Fraction of each query contig covered by minimap2 alignments to the target."""
    paf = Path(tmp) / "out.paf"
    with open(paf, "w") as out:
        subprocess.run([TOOLS / "minimap2", "-x", preset, "-t", "8", str(target), str(query_fasta)],
                       stdout=out, stderr=subprocess.DEVNULL, check=True)
    covered, length = defaultdict(int), {}
    for line in open(paf):
        f = line.split("\t")
        covered[f[0]] += int(f[3]) - int(f[2])
        length[f[0]] = int(f[1])
    return {q: min(1.0, covered[q] / length[q]) for q in length}


def coding_density(fasta, tmp):
    gff = Path(tmp) / "genes.gff"
    subprocess.run([TOOLS / "prodigal", "-p", "meta", "-g", "11", "-q", "-f", "gff", "-i", str(fasta),
                    "-o", str(gff)], check=True)
    coding = defaultdict(int)
    for line in open(gff):
        f = line.split("\t")
        if not line.startswith("#") and len(f) > 4:
            coding[f[0]] += int(f[4]) - int(f[3]) + 1
    return coding


def idxstats(bam):
    result = subprocess.run([TOOLS / "samtools", "idxstats", str(bam)], capture_output=True, text=True, check=True)
    return {f[0]: int(f[2]) for f in (line.split("\t") for line in result.stdout.splitlines())}


OUT.mkdir(parents=True, exist_ok=True)
meta = {f"{r['study']}__{r['library_id'].split(':')[-1]}": r["species_current"]
        for r in rows(ROOT / "manifest.csv", delimiter=",")}
grades = {v: {(r["sample_id"], r["target_id"]): r["grade"]
              for r in rows(RESULTS / "phase6_analysis" / v / "grades" / "bacterial_grades.tsv")}
          for v in ("v1", "v2")}
present_v1 = sorted(s for (s, t), g in grades["v1"].items() if t == MIXED and g in PRESENT)
orphans = {s for s in present_v1 if grades["v2"][(s, MERGED_INTO)] not in PRESENT}

with tempfile.TemporaryDirectory() as tmp:
    # Part 1: composition of the mixed bin.
    v1_fasta = RESULTS / "phase4_catalog" / "v1" / "bacterial_catalog.fna"
    mixed = read_fasta(v1_fasta, lambda n: n.startswith(MIXED + "|"))
    merged = read_fasta(v1_fasta, lambda n: n.startswith(MERGED_INTO + "|"))
    write_fasta(Path(tmp) / "mixed.fa", mixed)
    write_fasta(Path(tmp) / "merged.fa", merged)
    to_merged = aligned_fraction(Path(tmp) / "mixed.fa", Path(tmp) / "merged.fa", tmp)
    to_host = aligned_fraction(Path(tmp) / "mixed.fa", HOST, tmp)
    to_v2 = aligned_fraction(Path(tmp) / "mixed.fa", RESULTS / "phase4_catalog" / "v2" / "bacterial_catalog.fna", tmp)
    coding = coding_density(Path(tmp) / "mixed.fa", tmp)
    part = {c: "mycoplasma" if to_merged.get(c, 0) >= 0.5 else "other" for c in mixed}

    reads = defaultdict(lambda: defaultdict(int))  # library -> part -> reads
    contig_reads = defaultdict(int)
    for library in present_v1:
        counts = idxstats(RESULTS / "phase5_mapping" / "v1" / "bam" / "bacterial" / f"{library}.bam")
        for contig in mixed:
            reads[library][part[contig]] += counts.get(contig, 0)
            if library in orphans:
                contig_reads[contig] += counts.get(contig, 0)
    write("bac00033_contigs.tsv", [{
        "contig": c, "length": len(s), "gc": round((s.upper().count("G") + s.upper().count("C")) / len(s), 4),
        "coding_density": round(coding[c] / len(s), 4), "part": part[c],
        "aligned_to_bac00025": round(to_merged.get(c, 0), 3), "aligned_to_physalia": round(to_host.get(c, 0), 3),
        "aligned_to_v2_catalog": round(to_v2.get(c, 0), 3), "reads_in_orphan_libraries": contig_reads[c],
    } for c, s in mixed.items()])
    write("bac00033_libraries.tsv", [{
        "sample_id": s, "host": meta.get(s, ""), "orphan": s in orphans,
        "v2_bac00025_grade": grades["v2"][(s, MERGED_INTO)],
        "reads_mycoplasma_contigs": reads[s]["mycoplasma"], "reads_other_contigs": reads[s]["other"],
    } for s in present_v1])

    # Part 2: low-coding contigs in every v2 MAG representative.
    manifest = rows(RESULTS / "phase4_catalog" / "v2" / "bacterial_catalog.manifest.tsv")
    mags = {r["catalog_id"] for r in manifest if r["role"] == "mag" and r["status"] == "representative"}
    v2_fasta = RESULTS / "phase4_catalog" / "v2" / "bacterial_catalog.fna"
    contigs = read_fasta(v2_fasta, lambda n: n.split("|")[0] in mags)
    write_fasta(Path(tmp) / "mags.fa", contigs)
    mag_coding = coding_density(Path(tmp) / "mags.fa", tmp)
    low = {c for c, s in contigs.items() if mag_coding[c] / len(s) < 0.5}
    write_fasta(Path(tmp) / "low.fa", {c: contigs[c] for c in low})
    low_host = aligned_fraction(Path(tmp) / "low.fa", HOST, tmp) if low else {}
    mag_reads = defaultdict(lambda: defaultdict(int))
    for bam in sorted((RESULTS / "phase5_mapping" / "v2" / "bam" / "bacterial").glob("*.bam")):
        for contig, n in idxstats(bam).items():
            if contig in contigs:
                mag_reads[contig.split("|")[0]]["low" if contig in low else "coding"] += n
    screen = []
    for mag in sorted(mags):
        members = [c for c in contigs if c.split("|")[0] == mag]
        low_members = [c for c in members if c in low]
        total_reads = mag_reads[mag]["low"] + mag_reads[mag]["coding"]
        screen.append({
            "catalog_id": mag, "contigs": len(members), "length": sum(len(contigs[c]) for c in members),
            "low_coding_contigs": len(low_members), "low_coding_bp": sum(len(contigs[c]) for c in low_members),
            "low_coding_fraction": round(sum(len(contigs[c]) for c in low_members)
                                         / sum(len(contigs[c]) for c in members), 4),
            "low_coding_aligned_to_physalia": sum(low_host.get(c, 0) >= 0.5 for c in low_members),
            "reads_low_coding": mag_reads[mag]["low"],
            "read_fraction_low_coding": round(mag_reads[mag]["low"] / total_reads, 4) if total_reads else 0,
        })
    write("v2_mag_low_coding.tsv", screen)

# Summary of the quoted numbers.
contig_table = rows(OUT / "bac00033_contigs.tsv")
library_table = rows(OUT / "bac00033_libraries.tsv")
other = [r for r in contig_table if r["part"] == "other"]
myco = [r for r in contig_table if r["part"] == "mycoplasma"]


def density(records):
    return sum(float(r["coding_density"]) * int(r["length"]) for r in records) / sum(int(r["length"]) for r in records)


orphan_rows = [r for r in library_table if r["orphan"] == "True"]
lines = [
    f"BAC00033 contigs: {len(contig_table)}; aligned >=50% to BAC00025: {len(myco)}; other: {len(other)}",
    f"coding density: mycoplasma part {density(myco):.2f}, other part {density(other):.2f}",
    f"other contigs aligned >=50% to P. physalis: {sum(float(r['aligned_to_physalia']) >= 0.5 for r in other)}; "
    f"to v2 catalog: {sum(float(r['aligned_to_v2_catalog']) >= 0.5 for r in other)}",
    f"v1 libraries with BAC00033 present: {len(library_table)}; orphan libraries: {len(orphan_rows)}",
    f"orphan libraries: reads on mycoplasma contigs {sum(int(r['reads_mycoplasma_contigs']) for r in orphan_rows)}, "
    f"on other contigs {sum(int(r['reads_other_contigs']) for r in orphan_rows)}",
    f"orphan libraries with zero reads on mycoplasma contigs: "
    f"{sum(int(r['reads_mycoplasma_contigs']) == 0 for r in orphan_rows)}",
    f"other contigs with reads in orphan libraries: {sum(int(r['reads_in_orphan_libraries']) > 0 for r in other)}",
]
(OUT / "summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
for r in rows(OUT / "v2_mag_low_coding.tsv"):
    if int(r["low_coding_contigs"]):
        print("\t".join(r.values()))
