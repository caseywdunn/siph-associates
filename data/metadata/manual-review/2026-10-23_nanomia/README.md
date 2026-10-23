# Nanomia identifications (2026-10-23)

Six provisional `Nanomia sp.` identifications were resolved by matching museum
vouchers in Ahuja et al. (2026) to the sample names in their S8 Table. The published
S8 CO1 tree identifies NA28, WS3 and NA38 as *N. bijuga*, WS4 as *N. cara*, and NA29
as *Nanomia* sp. 1. Fig 4A identifies WS2 as *N. bijuga*. The earlier author-confirmed
NA37 identification remains unchanged.

- [Specimen identifications and evidence](identifications.tsv)
- [Published S8 tree](https://doi.org/10.1371/journal.pone.0351247.s016)
- [Published Fig 4 mitochondrial tree](https://doi.org/10.1371/journal.pone.0351247.g004)
- [Published S8 specimen table](https://doi.org/10.1371/journal.pone.0351247.s008)

The two downloaded figure images are retained alongside their SHA256 checksums
in `sources.json`. `before.json` records the input and reporting hashes before
this correction. Frozen sequence-analysis manifests and library identifiers
remain unchanged. The current determinations are applied through
`config/publication_host_names.tsv` and the documented publication metadata
overlay. This project correction does not update YPM or other public specimen
records.

The corrected 22-specimen Nanomia sample comprises 10 *N. bijuga*, 10 *N. septata*,
one *N. cara* and one *Nanomia* sp. 1. Four *N. bijuga* specimens have accepted
bacterial detections: NA31, NA33, NA37 and NA38. NA38 was previously pooled under
`Nanomia sp.`. All three original bacterial, viral and eukaryote detection tables
are byte-identical after rebuilding.

Regeneration uses the existing reporting-only workflows, in this order:

```bash
export XDG_CACHE_HOME=/tmp/siph-snakemake-cache
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
snakemake -s workflow/coauthor_reporting.smk --cores 1 --dry-run
snakemake -s workflow/coauthor_reporting.smk --cores 1 --printshellcmds
snakemake -s workflow/coauthor_metadata.smk all --cores 1 --dry-run
snakemake -s workflow/coauthor_metadata.smk all --cores 1 --printshellcmds
snakemake -s workflow/human_cutibacterium.smk --cores 1 --dry-run
snakemake -s workflow/human_cutibacterium.smk --cores 1 --printshellcmds
snakemake -s workflow/cutibacterium_sensitivity.smk --cores 1 --dry-run
snakemake -s workflow/cutibacterium_sensitivity.smk --cores 1 --printshellcmds
```

The executable is `/gpfs/gibbs/project/dunn/cwd7/conda_envs/snakemake/bin/snakemake`.
No raw-read processing, assembly, or mapping is scheduled by these entry points.
Flowcell tests are regenerated because host identification defines their strata;
community checks are refreshed against the corrected metadata. The manuscript
receives reviewed outputs through its checksum-validating asset sync.

Validation passed for all 213 libraries, the six exact voucher assignments,
the 10/10/1/1 Nanomia counts, and the four bacterial-positive N. bijuga
specimens. Frozen input hashes, the three detection tables, the community
and dispersion tests, and the Cutibacterium-omission community tests are
unchanged. See `verification.json`. The 35 relevant metadata tests pass.
Ruff lint passes; its format check reports two pre-existing formatting
differences in `summarize_dalila_revision.py`, which were left unchanged.
