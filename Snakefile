"""Publication workflow for the siphonophore-associate analysis."""
import csv
from pathlib import Path

configfile: "config/config.yaml"

ROOT = Path(workflow.basedir).resolve()
WORK = config["work_root"].rstrip("/")
SCRATCH = config["scratch_root"].rstrip("/")
FIXTURE_ROOT = config["fixture_root"].rstrip("/")

with open(config["samples"], newline="") as handle:
    SAMPLE_ROWS = list(csv.DictReader(handle, delimiter="\t"))
SAMPLES = {
    row["sample_id"]: row for row in SAMPLE_ROWS
    if row["include_primary"].lower() == "true"
}
SAMPLE_IDS = sorted(SAMPLES)

if not SAMPLE_IDS:
    raise WorkflowError("sample table has no included samples")
if len(SAMPLE_IDS) != len(set(SAMPLE_IDS)):
    raise WorkflowError("sample IDs are not unique")

wildcard_constraints:
    sample="[A-Za-z0-9_.-]+",
    route="P_physalis|N_septata",


include: "workflow/rules/common.smk"
include: "workflow/rules/preprocess.smk"
include: "workflow/rules/screens.smk"
include: "workflow/rules/host.smk"
include: "workflow/rules/validation.smk"
include: "workflow/rules/phase2.smk"
include: "workflow/rules/phase3.smk"
include: "workflow/rules/phase3_cohort.smk"
include: "workflow/rules/phase3_assembly.smk"
include: "workflow/rules/phase3_catalog.smk"
include: "workflow/rules/phase4_catalog.smk"
include: "workflow/rules/phase5_mapping.smk"
include: "workflow/rules/phase6_analysis.smk"
include: "workflow/rules/eukaryote_gate.smk"
include: "workflow/rules/mycoplasmatales_phylogeny.smk"


rule all:
    input:
        DEFAULT_TARGETS
