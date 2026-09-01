import os
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, atomic_move, ensure_parents, executable, sha256, version

ensure_parents(list(snakemake.output) + [snakemake.log[0]])
prefix = Path(str(snakemake.params.prefix))
bwa = executable(prefix, "bwa")
samtools = executable(prefix, "samtools")
jgi = executable(prefix, "jgi_summarize_bam_contig_depths")
metabat = executable(prefix, "metabat2")
environment = os.environ.copy()
environment["PATH"] = str(prefix / "bin") + os.pathsep + environment.get("PATH", "")
scratch_parent = Path(snakemake.config["scratch_root"]) / "phase3_benchmark" / "binning_tmp"
scratch_parent.mkdir(parents=True, exist_ok=True)

with tempfile.TemporaryDirectory(prefix=f"bin.{snakemake.wildcards.benchmark_id}.",
                                 dir=str(scratch_parent)) as temporary:
    temporary = Path(temporary)
    contigs = temporary / "assembly.fasta"
    os.link(snakemake.input.contigs, contigs)
    bam = temporary / "reads.sorted.bam"
    depth = temporary / "depth.tsv"
    bins_dir = temporary / "bins"
    bins_dir.mkdir()
    commands = []
    with open(snakemake.log[0], "w") as log:
        index_command = [bwa, "index", str(contigs)]
        commands.append(index_command)
        result = subprocess.run(index_command, stdout=log, stderr=subprocess.STDOUT, text=True, env=environment)
        if result.returncode:
            raise subprocess.CalledProcessError(result.returncode, index_command)
        mapping = [bwa, "mem", "-t", str(snakemake.threads), str(contigs),
                   str(snakemake.input.r1), str(snakemake.input.r2)]
        sorting = [samtools, "sort", "-@", str(max(1, snakemake.threads // 2)), "-o", str(bam), "-"]
        commands.extend([mapping, sorting])
        mapper = subprocess.Popen(mapping, stdout=subprocess.PIPE, stderr=log, env=environment)
        assert mapper.stdout is not None
        sorter = subprocess.Popen(sorting, stdin=mapper.stdout, stdout=log, stderr=log, env=environment)
        mapper.stdout.close()
        sorter_code = sorter.wait()
        mapper_code = mapper.wait()
        if mapper_code or sorter_code:
            raise RuntimeError(f"assembly backmap failed: bwa={mapper_code}, samtools={sorter_code}")
        depth_command = [jgi, "--outputDepth", str(depth), str(bam)]
        commands.append(depth_command)
        result = subprocess.run(depth_command, stdout=log, stderr=subprocess.STDOUT, text=True, env=environment)
        if result.returncode:
            raise subprocess.CalledProcessError(result.returncode, depth_command)
        bin_command = [metabat, "-i", str(contigs), "-a", str(depth), "-o", str(bins_dir / "bin"),
                       "-t", str(snakemake.threads), "-m", str(snakemake.params.min_contig),
                       "--seed", str(snakemake.params.seed)]
        commands.append(bin_command)
        result = subprocess.run(bin_command, stdout=log, stderr=subprocess.STDOUT, text=True, env=environment)
        if result.returncode:
            raise subprocess.CalledProcessError(result.returncode, bin_command)
    bins = sorted(bins_dir.glob("bin.*.fa"))
    archive = temporary / "bins.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        for path in bins:
            handle.add(path, arcname=path.name)
    total_bases = 0
    for path in bins:
        with path.open() as handle:
            total_bases += sum(len(line.strip()) for line in handle if not line.startswith(">"))
    atomic_move(depth, snakemake.output.depth)
    atomic_move(archive, snakemake.output.bins)

atomic_json(snakemake.output.metrics, {
    "benchmark_id": str(snakemake.wildcards.benchmark_id),
    "raw_bins": len(bins),
    "binned_bases": total_bases,
    "parameters": {"min_contig": int(snakemake.params.min_contig), "seed": int(snakemake.params.seed)},
    "software": {
        "bwa": version([bwa]), "samtools": version([samtools, "--version"]),
        "metabat2": version([metabat, "--help"]),
    },
    "assembly_provenance_sha256": sha256(snakemake.input.upstream),
})
