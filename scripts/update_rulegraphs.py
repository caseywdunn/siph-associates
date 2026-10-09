#!/usr/bin/env python3
"""Generate or check documentation graphs without executing workflow jobs."""

import argparse
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs" / "workflow"
SNAKEMAKE_VERSION = "7.24.0"
GRAPHVIZ_VERSION = "2.40.1"


def run(command: list[str], env: dict[str, str]) -> bytes:
    """Capture a command's output, keeping diagnostics separate from DOT."""
    return subprocess.run(
        command, cwd=ROOT, env=env, check=True, stdout=subprocess.PIPE
    ).stdout


def graph_command(spec: dict, snakemake: str, mode: str) -> list[str]:
    """Use the same entry point, target and configuration for graph and preview."""
    command = [
        snakemake,
        *spec["targets"],
        "--snakefile",
        "Snakefile",
        "--configfile",
        spec["config"],
        "--cores",
        "1",
        "--nolock",
        "--rerun-triggers",
        "mtime",
        mode,
    ]
    if spec.get("allowed_rules"):
        command.extend(["--allowed-rules", *spec["allowed_rules"]])
    return command


def main() -> None:
    """Render all requested graphs, or fail if checked-in graphs are stale."""
    specs = json.loads((DOCS / "graphs.json").read_text())
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("graphs", nargs="*", help="Graph IDs; default: all")
    parser.add_argument("--check", action="store_true", help="Compare without writing")
    parser.add_argument(
        "--require-tracked", action="store_true", help="Also require Git-tracked images"
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Preview these scopes; do not render"
    )
    args = parser.parse_args()
    unknown = set(args.graphs) - specs.keys()
    if unknown:
        parser.error(f"Unknown graph IDs: {', '.join(sorted(unknown))}")
    if args.dry_run and (args.check or args.require_tracked):
        parser.error("--dry-run cannot be combined with graph checks")
    selected = args.graphs or list(specs)
    snakemake = os.environ.get("SNAKEMAKE", "snakemake")
    dot = os.environ.get("DOT", "dot")
    env = {**os.environ, "PYTHONHASHSEED": "0"}
    env.setdefault("XDG_CACHE_HOME", str(ROOT / ".cache"))
    # A user's default profile must not silently change the documented scope.
    env.pop("SNAKEMAKE_PROFILE", None)
    version = run([snakemake, "--version"], env).decode().strip()
    if version != SNAKEMAKE_VERSION:
        raise SystemExit(f"Use Snakemake {SNAKEMAKE_VERSION}; found {version}")
    if not args.dry_run:
        version = subprocess.run(
            [dot, "-V"], check=True, capture_output=True, text=True
        )
        if (
            f"graphviz version {GRAPHVIZ_VERSION} "
            not in version.stderr + version.stdout
        ):
            raise SystemExit(f"Use Graphviz {GRAPHVIZ_VERSION}")

    destination = DOCS / "graphs"
    # Stage every graph before replacing anything: a failed Snakemake or Graphviz
    # call must not damage the previous documentation.
    with tempfile.TemporaryDirectory(prefix="siph-rulegraphs-") as temporary:
        staging = Path(temporary)
        expected = []
        for name in selected:
            spec = specs[name]
            print(f"{name}: {' '.join(spec['targets'])}", flush=True)
            mode = "--dry-run" if args.dry_run else "--rulegraph"
            data = run(graph_command(spec, snakemake, mode), env)
            if args.dry_run:
                print(data.decode(), end="", flush=True)
                continue
            (staging / f"{name}.dot").write_bytes(data)
            run(
                [
                    dot,
                    "-Tsvg",
                    str(staging / f"{name}.dot"),
                    "-o",
                    str(staging / f"{name}.svg"),
                ],
                env,
            )
            expected.extend([f"{name}.dot", f"{name}.svg"])
        if args.require_tracked:
            run(
                ["git", "ls-files", "--error-unmatch", "--"]
                + [str((destination / name).relative_to(ROOT)) for name in expected],
                env,
            )
        if args.check:
            stale = [
                name
                for name in expected
                if not (destination / name).is_file()
                or (destination / name).read_bytes() != (staging / name).read_bytes()
            ]
            if stale:
                raise SystemExit("Missing or stale rule graphs: " + ", ".join(stale))
            print(f"All {len(selected)} rule graphs match.")
        elif not args.dry_run:
            destination.mkdir(parents=True, exist_ok=True)
            for name in expected:
                # The temporary replacement lives beside its destination so rename
                # stays atomic even when /tmp and the checkout are on different disks.
                with tempfile.NamedTemporaryFile(
                    dir=destination, delete=False
                ) as handle:
                    replacement = Path(handle.name)
                try:
                    shutil.copyfile(staging / name, replacement)
                    replacement.chmod(0o644)
                    replacement.replace(destination / name)
                finally:
                    replacement.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
