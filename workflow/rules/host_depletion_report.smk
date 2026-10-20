"""Small mapping-summary inputs only; no FASTQ access or upstream construction."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path("scripts").resolve()))
from summarize_host_depletion import SOURCE_FILES, report_inputs

HOST_DEPLETION_REPORT_CONFIG = "config/host_depletion_report.json"
HOST_DEPLETION_REPORT = json.loads(Path(HOST_DEPLETION_REPORT_CONFIG).read_text())
HOST_DEPLETION_REPORT_OUT = HOST_DEPLETION_REPORT["output_dir"]
HOST_DEPLETION_REQUIRE_COMPLETE = str(
    config.get("host_depletion_require_complete", False)
).lower() == "true"


rule report_host_depletion:
    input:
        configuration=HOST_DEPLETION_REPORT_CONFIG,
        sources=SOURCE_FILES,
        records=lambda w: report_inputs(
            HOST_DEPLETION_REPORT, HOST_DEPLETION_REQUIRE_COMPLETE
        ),
    output:
        table=f"{HOST_DEPLETION_REPORT_OUT}/mapping_stats.tsv",
        markdown=f"{HOST_DEPLETION_REPORT_OUT}/mapping_stats.md",
        provenance=f"{HOST_DEPLETION_REPORT_OUT}/provenance.json",
        source_snapshot=[
            f"{HOST_DEPLETION_REPORT_OUT}/provenance/source/{p}"
            for p in [*SOURCE_FILES, HOST_DEPLETION_REPORT_CONFIG]
        ],
    params:
        python=sys.executable,
        require_complete="--require-complete" if HOST_DEPLETION_REQUIRE_COMPLETE else "",
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
    shell:
        """
        {params.python:q} scripts/summarize_host_depletion.py \
          --config {input.configuration:q} {params.require_complete}
        """
