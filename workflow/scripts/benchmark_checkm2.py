import csv
import os
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, atomic_move, ensure_parents, executable, sha256, version
from checkm2_result import classify_checkm2_exit

ensure_parents(list(snakemake.output) + [snakemake.log[0]])
prefix = Path(str(snakemake.params.prefix))
checkm2 = executable(prefix, "checkm2")
environment = os.environ.copy()
environment["PATH"] = str(prefix / "bin") + os.pathsep + environment.get("PATH", "")

with tempfile.TemporaryDirectory(prefix=f"checkm2.{snakemake.wildcards.benchmark_id}.",
                                 dir=str(Path(snakemake.output.quality).parent)) as temporary:
    temporary = Path(temporary)
    bins_dir, output_dir = temporary / "bins", temporary / "output"
    bins_dir.mkdir()
    with tarfile.open(snakemake.input.bins, "r:gz") as archive:
        members = archive.getmembers()
        if any(member.name.startswith("/") or ".." in Path(member.name).parts for member in members):
            raise ValueError("unsafe path in bin archive")
        archive.extractall(bins_dir, members=members)
    bins = sorted(bins_dir.glob("*.fa"))
    quality = temporary / "quality_report.tsv"
    command = None
    status = "no_bins"
    with open(snakemake.log[0], "w") as log:
        if bins:
            command = [checkm2, "predict", "--threads", str(snakemake.threads), "--input", str(bins_dir),
                       "--extension", "fa", "--database_path", str(snakemake.input.database),
                       "--output-directory", str(output_dir), "--force"]
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True, env=environment)
            log.flush()
            status = classify_checkm2_exit(result.returncode, Path(snakemake.log[0]).read_text())
            if status == "assessed":
                produced = output_dir / "quality_report.tsv"
                if not produced.is_file():
                    raise FileNotFoundError("CheckM2 did not produce quality_report.tsv")
                os.replace(produced, quality)
            else:
                quality.write_text("Name\tCompleteness\tContamination\n")
                log.write("Recorded bins as unassessable because CheckM2 found no DIAMOND annotations.\n")
        else:
            quality.write_text("Name\tCompleteness\tContamination\n")
            log.write("No MetaBAT2 bins; CheckM2 not run.\n")
    rows = list(csv.DictReader(quality.open(), delimiter="\t"))
    high = sum(float(row["Completeness"]) >= 90 and float(row["Contamination"]) < 5 for row in rows)
    medium = sum(float(row["Completeness"]) >= 50 and float(row["Contamination"]) < 10 for row in rows)
    atomic_move(quality, snakemake.output.quality)

atomic_json(snakemake.output.metrics, {
    "benchmark_id": str(snakemake.wildcards.benchmark_id),
    "checkm2_status": status,
    "assessed_bins": len(rows),
    "unassessed_bins": len(bins) - len(rows),
    "medium_quality_or_better_bins": medium,
    "high_quality_bins": high,
    "command": command,
    "software": version([checkm2, "--version"]),
    "bin_metrics_sha256": sha256(snakemake.input.upstream),
})
