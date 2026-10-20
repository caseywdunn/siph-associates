"""Report existing host mappings without permitting sequence analyses to rerun."""

include: "rules/host_depletion_report.smk"


rule all:
    input: rules.report_host_depletion.output
    default_target: True
