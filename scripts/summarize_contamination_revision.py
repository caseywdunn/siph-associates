#!/usr/bin/env python3
"""Report both exclusion tiers, retaining the full cohort in the contamination screen."""
from __future__ import annotations
import csv
import hashlib
import json
import os
from collections import Counter, defaultdict
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
from contamination_policy import ContaminationPolicy, POLICY_PATHS, read_tsv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = Path("figures/contamination")
out.mkdir(parents=True, exist_ok=True)
settings = json.loads(Path("config/reporting_sources.json").read_text())
inputs = [Path("config/reporting_sources.json"), *map(Path, POLICY_PATHS),
          Path("workflow/scripts/contamination_policy.py"), Path(__file__),
          Path("figures/human_screen/sample_summary.tsv"),
          Path(settings["analysis_root"]) / "grades/bacterial_grades.tsv",
          Path("figures/cohort/eukaryote_evidence.tsv"),
          Path("figures/cohort/sample_metadata.tsv"), Path(settings["catalog_manifest"])]
catalog = read_tsv(inputs[-1])
policy = ContaminationPolicy(catalog=catalog)
bacteria = [policy.annotate(r) for r in read_tsv(inputs[6])]
eukaryotes = read_tsv(inputs[7])
metadata = {r["analysis_sample_id"]: r for r in read_tsv(inputs[8])}
human = read_tsv(inputs[5])
acnes = {r["sample_id"]: r for r in bacteria if "s__Cutibacterium acnes" in r["taxonomy"]}
assert len(human) == len(metadata) == 213
rows = []
for r in human:
    sample = r["sample_id"]
    bacterial = acnes[sample]
    rows.append({"sample_id": sample, "host": metadata[sample]["species_current"],
        "screened_pairs": int(r["raw_selected_pairs"]),
        "human_compatible_pairs": int(r["stringent_human_compatible_pairs"]),
        "human_pairs_per_million": 1e6 * int(r["stringent_human_compatible_pairs"]) / int(r["raw_selected_pairs"]),
        "cutibacterium_acnes_breadth": float(bacterial["breadth"]),
        "excluded_library": sample in policy.excluded,
        "reason": policy.excluded.get(sample, {}).get("reason", "")})
def write(name, records, fields=None):
    with (out / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
write("sample_screen.tsv", rows)
plt.rcParams.update({"pdf.fonttype": 42, "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(6.4, 4.8))
for excluded, color, label in [(False, "#0072B2", "Retained specimens"), (True, "#D55E00", "Excluded specimens")]:
    subset = [r for r in rows if r["excluded_library"] == excluded]
    ax.scatter([r["human_pairs_per_million"] for r in subset], [100*r["cutibacterium_acnes_breadth"] for r in subset], s=22, color=color, alpha=.75, label=label)
    if excluded:
        for r in subset:
            label = r["sample_id"].split("__", 1)[1]
            placement = {"xytext": (5, 5)}
            if label == "D666-D10":
                placement = {"xytext": (-5, -5), "ha": "right", "va": "top"}
            elif label == "NA31":
                placement = {"xytext": (5, 5), "ha": "left", "va": "bottom"}
            ax.annotate(label, (r["human_pairs_per_million"], 100*r["cutibacterium_acnes_breadth"]),
                        textcoords="offset points", fontsize=8, **placement)
ax.axvline(500, ls="--", lw=.8, color="gray")
ax.axhline(50, ls="--", lw=.8, color="gray")
ax.set_xscale("symlog", linthresh=2)
ax.set_xlim(0, max(r["human_pairs_per_million"] for r in rows) * 1.7)
ax.set_xlabel("Human-compatible read pairs per million screened pairs")
ax.set_ylabel(r"$\mathit{Cutibacterium\ acnes}$ genome breadth (%)")
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(out / "human_cutibacterium.pdf", bbox_inches="tight")
fig.savefig(out / "human_cutibacterium.png", bbox_inches="tight", dpi=180)
plt.close(fig)
accepted = [r for r in bacteria if r["role"] != "decoy" and r["grade"] in {"validated", "high_confidence"}]
contaminants = [r for r in accepted if r["contaminant_taxon"] == "true"]
write("bacterial_contaminant_detections.tsv", contaminants)
ce = [r for r in eukaryotes if r["contaminant_taxon"] == "true"]
write("eukaryote_contaminant_detections.tsv", ce)
counts = Counter((r["target_id"], r["taxonomy"], r["contaminant_group"]) for r in contaminants)
write("contaminant_targets.tsv", [dict(target_id=key[0],taxonomy=key[1],category=key[2],detections=n) for key,n in sorted(counts.items())])
summary = {"screened_libraries": len(rows), "excluded_libraries": len(policy.excluded),
    "retained_libraries": len(rows)-len(policy.excluded),
    "excluded_ids": sorted(policy.excluded), "bacterial_contaminant_detections": len(contaminants),
    "bacterial_contaminant_positive_libraries": len({r["sample_id"] for r in contaminants}),
    "eukaryote_contaminant_detections": len(ce),
    "eukaryote_contaminant_quality": "Both libraries were scored moderate quality by Church et al. (2025)",
    "decision_basis": "Author thresholds chosen after inspection at a gap in the distributions; no extraction blanks.",
    "grade_codes": "The stable internal high_confidence code is displayed as assembly-supported (bacteria/viruses) or reads + assembly (eukaryotes)."}
(out / "summary.json").write_text(json.dumps(summary,indent=2)+"\n")
outputs = sorted(p for p in out.iterdir() if p.is_file() and p.name != "provenance.json")
record = {"status":"complete", "run_id":os.environ.get("SLURM_JOB_ID", "unknown"),
    "inputs_sha256": {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    "outputs_sha256": {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs}}
(out / "provenance.json").write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps(summary,indent=2))
