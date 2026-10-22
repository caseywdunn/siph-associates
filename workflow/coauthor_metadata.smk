"""Tissue/depth summaries and conditional community tests from completed evidence.

This bounded entry point has no sequence-processing rules. Missing mapping or
evidence tables fail rather than scheduling a read-processing round.
"""

include: "rules/coauthor_metadata.smk"


rule all:
    input: rules.record_coauthor_metadata.output
    default_target: True
