"""Interpret CheckM2 exits without masking general execution failures."""

NO_ANNOTATIONS_MESSAGE = "ERROR: No DIAMOND annotation was generated. Exiting"


def classify_checkm2_exit(returncode, log_text):
    if returncode == 0:
        return "assessed"
    if NO_ANNOTATIONS_MESSAGE in log_text:
        return "no_annotations"
    raise RuntimeError(f"CheckM2 failed with exit code {returncode}")
