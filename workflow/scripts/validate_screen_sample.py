import csv
import json
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, ensure_parents, fastq_pair_count, file_record

ensure_parents([snakemake.output[0], snakemake.log[0]])
errors = []
sample = str(snakemake.wildcards.sample)


def read_header(path, delimiter="\t"):
    with open(path, newline="") as handle:
        return next(csv.reader(handle, delimiter=delimiter), [])


try:
    trimmed_pairs = fastq_pair_count(snakemake.input.trim_r1, snakemake.input.trim_r2)
except Exception as exc:
    errors.append(f"trimmed FASTQ pair validation: {exc}")
    trimmed_pairs = None

try:
    fastp = json.loads(Path(snakemake.input.fastp).read_text())
    fastp_pairs = int(fastp["summary"]["after_filtering"]["total_reads"]) // 2
    if int(fastp["summary"]["after_filtering"]["total_reads"]) % 2:
        errors.append("fastp after-filtering read count is odd")
    if trimmed_pairs is not None and fastp_pairs != trimmed_pairs:
        errors.append(f"fastp/output pair mismatch: report={fastp_pairs}, files={trimmed_pairs}")
except Exception as exc:
    errors.append(f"fastp report validation: {exc}")
    fastp_pairs = None

try:
    trim_provenance = json.loads(Path(snakemake.input.trim_provenance).read_text())
    expected_entering = min(int(snakemake.params.raw_pairs), int(snakemake.params.cap_pairs))
    if trim_provenance.get("sample_id") != sample:
        errors.append("trim provenance sample ID mismatch")
    if int(trim_provenance["raw_pairs_expected"]) != int(snakemake.params.raw_pairs):
        errors.append("trim provenance raw-pair count mismatch")
    if int(trim_provenance["pairs_entering_fastp"]) != expected_entering:
        errors.append("trim provenance cap count mismatch")
    if fastp_pairs is not None and int(trim_provenance["pairs_after_fastp"]) != fastp_pairs:
        errors.append("trim provenance post-fastp count mismatch")
except Exception as exc:
    errors.append(f"trim provenance validation: {exc}")

expected_headers = {
    "genus": ["name", "taxonomy_id", "taxonomy_lvl", "kraken_assigned_reads", "added_reads", "new_est_reads", "fraction_total_reads"],
    "species": ["name", "taxonomy_id", "taxonomy_lvl", "kraken_assigned_reads", "added_reads", "new_est_reads", "fraction_total_reads"],
    "sylph": ["Sample_file", "Genome_file", "Taxonomic_abundance", "Sequence_abundance", "Adjusted_ANI", "True_cov", "ANI_5-95_percentile", "Eff_lambda", "Lambda_5-95_percentile", "Median_cov", "Mean_cov_geq1", "Containment_ind", "Naive_ANI", "kmers_reassigned", "Contig_name"],
}
row_counts = {}
for label, path in (("genus", snakemake.input.genus), ("species", snakemake.input.species), ("sylph", snakemake.input.sylph)):
    try:
        header = read_header(path)
        if header != expected_headers[label]:
            errors.append(f"{label} header mismatch: {header}")
        with open(path, newline="") as handle:
            row_counts[label] = sum(1 for _ in csv.DictReader(handle, delimiter="\t"))
    except Exception as exc:
        errors.append(f"{label} table validation: {exc}")

try:
    kraken_rows = 0
    with open(snakemake.input.kraken) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 8:
                errors.append(f"Kraken2 report row has {len(fields)} columns, expected 8")
                break
            kraken_rows += 1
    if kraken_rows == 0:
        errors.append("Kraken2 report has no rows")
    row_counts["kraken"] = kraken_rows
except Exception as exc:
    errors.append(f"Kraken2 report validation: {exc}")

try:
    with tarfile.open(snakemake.input.phyloflash, "r:gz") as outer:
        names = [member.name for member in outer.getmembers() if member.isfile()]
        nested = [name for name in names if name.endswith(".phyloFlash.tar.gz")]
        if not nested:
            errors.append("phyloFlash archive lacks nested result archive")
        if not any(name.endswith(".phyloFlash.html") for name in names):
            errors.append("phyloFlash archive lacks HTML report")
        if nested:
            stream = outer.extractfile(nested[0])
            if stream is None:
                errors.append("cannot read nested phyloFlash result archive")
            else:
                with tarfile.open(fileobj=stream, mode="r:gz") as inner:
                    inner_names = [member.name for member in inner.getmembers() if member.isfile()]
                if not any(name.endswith(".phyloFlash.report.csv") for name in inner_names):
                    errors.append("phyloFlash nested archive lacks run report")
                if not any(name.endswith(".phyloFlash.NTUfull_abundance.csv") for name in inner_names):
                    errors.append("phyloFlash nested archive lacks full NTU abundance")
    row_counts["phyloflash_outer_files"] = len(names)
except Exception as exc:
    errors.append(f"phyloFlash archive validation: {exc}")

for label, path in (("Kraken2", snakemake.input.kraken_provenance),
                    ("sylph", snakemake.input.sylph_provenance),
                    ("phyloFlash", snakemake.input.phyloflash_provenance)):
    try:
        record = json.loads(Path(path).read_text())
        if record.get("sample_id") != sample:
            errors.append(f"{label} provenance sample ID mismatch")
    except Exception as exc:
        errors.append(f"{label} provenance validation: {exc}")

Path(snakemake.log[0]).write_text("\n".join([f"sample={sample}", f"errors={len(errors)}", *errors]) + "\n")
if errors:
    raise ValueError("; ".join(errors))

atomic_json(snakemake.output[0], {
    "status": "PASS",
    "sample_id": sample,
    "library_id": str(snakemake.params.library_id),
    "specimen_id": str(snakemake.params.specimen_id),
    "study": str(snakemake.params.study),
    "host_route": str(snakemake.params.route),
    "raw_pairs": int(snakemake.params.raw_pairs),
    "cap_pairs": int(snakemake.params.cap_pairs),
    "cap_applied": int(snakemake.params.raw_pairs) > int(snakemake.params.cap_pairs),
    "pairs_entering_fastp": min(int(snakemake.params.raw_pairs), int(snakemake.params.cap_pairs)),
    "pairs_after_fastp": trimmed_pairs,
    "screen_rows": row_counts,
    "validated_outputs": [file_record(path) for path in list(snakemake.input)],
})
