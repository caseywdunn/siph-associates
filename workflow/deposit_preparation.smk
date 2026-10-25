"""Bounded deposit preparation from reviewed catalog FASTAs; no read processing."""
import csv
import json
from pathlib import Path

include: "rules/deposit_preparation.smk"

rule prepare_deposit:
    input: rules.deposit_preparation.input
    default_target: True
