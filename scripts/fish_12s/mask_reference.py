#!/usr/bin/env python3
"""Make a DUST-masked search derivative without changing the source references."""

import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO

from pilot import digest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dustmasker", required=True)
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=False)
    soft = out / "softmasked.fasta"
    command = [
        args.dustmasker,
        "-in",
        str(args.input),
        "-out",
        str(soft),
        "-outfmt",
        "fasta",
        "-window",
        "64",
        "-level",
        "20",
        "-linker",
        "1",
    ]
    subprocess.run(command, check=True)
    source = {r.id: str(r.seq) for r in SeqIO.parse(args.input, "fasta")}
    masked_bases, total_bases, seen = 0, 0, set()
    with (out / "fish_12s_masked.fasta").open("w") as handle:
        for record in SeqIO.parse(soft, "fasta"):
            seq = str(record.seq)
            if record.id in seen or seq.upper() != source[record.id]:
                raise ValueError(f"DUST changed reference identity: {record.id}")
            seen.add(record.id)
            masked_bases += sum(base.islower() for base in seq)
            total_bases += len(seq)
            handle.write(
                f">{record.id}\n{''.join('N' if b.islower() else b for b in seq)}\n"
            )
    if seen != set(source):
        raise ValueError("DUST dropped reference records")
    run = {
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "command": command,
        "script_sha256": digest(Path(__file__)),
        "input_sha256": digest(args.input),
        "dustmasker_version": subprocess.check_output(
            [args.dustmasker, "-version"], text=True
        ).strip(),
        "sequences": len(seen),
        "total_bases": total_bases,
        "masked_bases": masked_bases,
        "outputs": {str(p): digest(p) for p in out.iterdir() if p.is_file()},
    }
    (out / "run.json").write_text(json.dumps(run, indent=2) + "\n")
    print(json.dumps({k: run[k] for k in ["sequences", "total_bases", "masked_bases"]}))


if __name__ == "__main__":
    main()
