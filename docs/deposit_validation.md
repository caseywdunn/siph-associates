# Local checks for deposit sequences

The validation launchers keep software, temporary files and results on GPFS.
They neither submit sequences nor supply author contact information.

Install table2asn and the ribotyper component of Ribovore from their official
sources with `sbatch -p ycga scripts/install_deposit_validation_tools.sbatch`.
The isolated installation is `data/software/deposit_validation/`; version/help
output, downloaded installer and checksums are recorded per execution. The
Ribovore installer is checked for release 1.0.5 and tested on its official example.
No existing conda environment or shell configuration is changed.

After the candidate bundle is prepared, run
`sbatch -p ycga scripts/deposit_sequence_validation.sbatch`. It checks the 65
packaged non-host SSUs and corrected MAG 16S queries with ribotyper. It also asks
table2asn to generate validation and discrepancy reports for the candidate MAGs.
If table2asn cannot produce reports without a real submission template, the
launcher stops and retains its diagnostic log. If reports are produced, they
remain structural drafts: a missing author template is recorded explicitly and
no output is labeled submission-ready. Any FCS sequence change requires a new
validation execution on the final bundle.

The table2asn `.sbt` submission template is genuinely absent. It must contain
real author contact and sequencing-credit information. Viral, SSU and marker
source and feature tables also remain incomplete; their required annotations
cannot be inferred from the FASTA identifiers. PGAP installation and taxcheck
are separate pending work.

Sources: [NCBI table2asn documentation](https://www.ncbi.nlm.nih.gov/genbank/table2asn/),
[NCBI table2asn README](https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/DOCUMENTATION/table2asn_readme.txt),
[Ribovore installation](https://github.com/NLM-DIR/ribovore/blob/master/documentation/install.md),
[ribotyper usage](https://github.com/NLM-DIR/ribovore/blob/master/documentation/ribotyper.md).
