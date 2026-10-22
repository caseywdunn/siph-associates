"""Count stringent human-compatible read pairs in query-grouped alignments."""

import argparse
from collections import Counter
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path

import pysam


def passes_read(read, minimum_length, minimum_identity, minimum_mapq):
    """Require a unique, near-exact, full-length primary mate alignment."""
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    if read.mapping_quality < minimum_mapq or read.query_length < minimum_length:
        return False
    if not read.has_tag("NM") or not read.cigartuples:
        return False
    # Clipping is prohibited even if an alternative aligner emits local matches.
    if any(op in (4, 5) for op, _ in read.cigartuples):
        return False
    aligned_columns = sum(
        length for op, length in read.cigartuples if op in (0, 1, 2, 7, 8)
    )
    if (
        not aligned_columns
        or 1 - read.get_tag("NM") / aligned_columns < minimum_identity
    ):
        return False
    seq = read.query_sequence.upper()
    proportions = [count / len(seq) for count in Counter(seq).values()]
    entropy = -sum(p * math.log2(p) for p in proportions)
    return entropy >= 1.2 and set(seq) <= set("ACGT")


def passes_pair(reads, minimum_length=100, minimum_identity=0.99, minimum_mapq=30):
    if len(reads) != 2:
        return False
    first, second = sorted(reads, key=lambda r: not r.is_read1)
    return (
        first.is_read1
        and second.is_read2
        and all(r.is_paired and r.is_proper_pair for r in reads)
        and first.reference_id == second.reference_id
        and first.is_reverse != second.is_reverse
        and 0 < abs(first.template_length) <= 1000
        and all(
            passes_read(r, minimum_length, minimum_identity, minimum_mapq)
            for r in reads
        )
    )


def assess(bam, config, sample, input_metrics, fastp_metrics, output):
    output.mkdir(parents=True, exist_ok=True)
    aligned_pairs, retained_pairs = 0, 0
    chromosomes, windows, sequence_pairs = set(), set(), set()
    with (
        pysam.AlignmentFile(bam, "rb") as alignments,
        gzip.open(output / "candidates_R1.fastq.gz", "wt") as r1,
        gzip.open(output / "candidates_R2.fastq.gz", "wt") as r2,
        open(output / "candidates.fasta", "w") as fasta,
        open(output / "pairs.tsv", "w") as table,
    ):
        table.write(
            "pair_id\tchromosome\tstart\tend\tmate1_nm\tmate2_nm\tmate1_length\tmate2_length\tsequence_sha256\n"
        )
        for name, group in itertools.groupby(alignments, key=lambda r: r.query_name):
            reads = [
                read
                for read in group
                if not read.is_secondary and not read.is_supplementary
            ]
            if any(not r.is_unmapped for r in reads):
                aligned_pairs += 1
            if not passes_pair(
                reads,
                config["minimum_length"],
                config["minimum_identity"],
                config["minimum_mapq"],
            ):
                continue
            retained_pairs += 1
            first, second = sorted(reads, key=lambda r: not r.is_read1)
            pair_id = f"{sample}__{name}"
            sequences = [first.get_forward_sequence(), second.get_forward_sequence()]
            digest = hashlib.sha256(
                (sequences[0] + ":" + sequences[1]).encode()
            ).hexdigest()
            sequence_pairs.add(digest)
            for mate, (read, handle, sequence) in enumerate(
                zip((first, second), (r1, r2), sequences, strict=True), 1
            ):
                quality = "".join(chr(q + 33) for q in read.get_forward_qualities())
                handle.write(f"@{pair_id}\n{sequence}\n+\n{quality}\n")
                fasta.write(f">{pair_id}/{mate}\n{sequence}\n")
            start = min(r.reference_start for r in reads)
            end = max(r.reference_end for r in reads)
            chromosome = first.reference_name
            chromosomes.add(chromosome)
            windows.add((chromosome, start // 1000000))
            table.write(
                f"{pair_id}\t{chromosome}\t{start}\t{end}\t{first.get_tag('NM')}\t{second.get_tag('NM')}\t{first.query_length}\t{second.query_length}\t{digest}\n"
            )
    before = input_metrics["selected_pairs"]
    after = fastp_metrics["summary"]["after_filtering"]["total_reads"] // 2
    summary = {
        "sample_id": sample,
        "selection": input_metrics["selection"],
        "raw_selected_pairs": before,
        "qc_pairs": after,
        "bowtie2_concordant_pairs": aligned_pairs,
        "stringent_human_compatible_pairs": retained_pairs,
        "distinct_sequence_pairs": len(sequence_pairs),
        "reference_sequences": len(chromosomes),
        "occupied_1Mb_windows": len(windows),
        "percent_selected_raw_pairs": 100 * retained_pairs / before if before else None,
        "percent_qc_pairs": 100 * retained_pairs / after if after else None,
        "taxonomic_confirmation": "pending; human-compatible mapping alone",
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("config", "sample", "bam", "input_metrics", "fastp_metrics", "output"):
        parser.add_argument("--" + name.replace("_", "-"), required=True)
    args = parser.parse_args()
    assess(
        args.bam,
        json.loads(Path(args.config).read_text()),
        args.sample,
        json.loads(Path(args.input_metrics).read_text()),
        json.loads(Path(args.fastp_metrics).read_text()),
        Path(args.output),
    )


if __name__ == "__main__":
    main()
