"""Record exact reference files and mapper indexes used for human screening."""

import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    paths = [config["human_fasta"], *config["negative_references"].values()]
    paths += [
        config["human_index"] + suffix
        for suffix in (
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        )
    ]
    records = []
    for name in paths:
        path = Path(name)
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(block)
        records.append(
            {
                "path": str(path.resolve()),
                "bytes": path.stat().st_size,
                "sha256": digest.hexdigest(),
            }
        )
    args.output.write_text(json.dumps(records, indent=2) + "\n")


if __name__ == "__main__":
    main()
