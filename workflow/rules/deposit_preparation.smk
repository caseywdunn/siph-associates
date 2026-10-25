"""
Curate the reviewed MAGs and stage sequence sets for public deposition.

1. Apply recorded contig removals, terminal adaptor trims and the palindrome
   repair: prepare_deposit_mag.
2. Reassess every edited MAG: assess_deposit_mag_checkm2.
3. Package all sequences, candidate/archival manifests and draft source metadata:
   package_deposit_sequences.
4. Screen every sequence set for adaptors: screen_deposit_adaptors.
5. Explicitly opt in to the large-memory cross-species screen:
   screen_deposit_mag_fcsgx, through target deposit_fcsgx.

Review preparation and adaptor reports before enabling FCS-GX. Review its reports before
final table2asn, PGAP taxcheck, ribovore and metadata validation; no submission is
made by this workflow. Original mapping catalog sequences remain unchanged.
"""

with open("config/deposit_preparation.json") as handle:
    DEPOSIT_CONFIG = json.load(handle)
DEP = DEPOSIT_CONFIG["output_root"]
DEP_PY = DEPOSIT_CONFIG["tools"]["python"]
with open(DEPOSIT_CONFIG["decisions"], newline="") as handle:
    DEPOSIT_MAGS = {row["mag_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
DEPOSIT_EDITED = set()
for table in ("removals", "trims", "repairs"):
    with open(DEPOSIT_CONFIG[table], newline="") as handle:
        DEPOSIT_EDITED.update(row["mag_id"] for row in csv.DictReader(handle, delimiter="\t"))
DEPOSIT_SET_IDS = ["mags", "viruses", "viral_endogenous", "ssu", "markers"]


def deposit_taxid(wildcards):
    with open(DEPOSIT_CONFIG["organism_names"], newline="") as handle:
        rows = {r["mag_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    taxid = rows[wildcards.mag]["screening_taxid"]
    if not taxid.isdigit():
        raise WorkflowError(f"No verified NCBI screening taxid for {wildcards.mag}")
    return taxid


def deposit_division_options(wildcards):
    taxid = deposit_taxid(wildcards)
    with open(DEPOSIT_CONFIG["fcs_gx_division_overrides"], newline="") as handle:
        overrides = {r["ncbi_taxid"]: r for r in csv.DictReader(handle, delimiter="\t")}
    if taxid not in overrides:
        return []
    row = overrides[taxid]
    if row["asserted_gx_division"] != "prok:mycoplasmas" or taxid != "2979199":
        raise WorkflowError("Unreviewed FCS-GX division override")
    return ["--div", row["asserted_gx_division"]]


rule prepare_deposit_mag:
    input:
        script="workflow/scripts/prepare_deposit_mag.py",
        fasta=lambda w: DEPOSIT_MAGS[w.mag]["source_fasta"],
        removals=DEPOSIT_CONFIG["removals"],
        trims=DEPOSIT_CONFIG["trims"],
        repairs=DEPOSIT_CONFIG["repairs"],
    output:
        fasta=f"{DEP}/cleaned_mags/{{mag}}.fa",
        report=f"{DEP}/edits/{{mag}}.json",
    wildcard_constraints:
        mag="|".join(DEPOSIT_MAGS),
    threads: 1
    resources: mem_mb=2000, runtime=15, partition="ycga"
    log: f"{DEP}/logs/prepare/{{mag}}.log"
    shell:
        """
        {DEP_PY:q} {input.script:q} --mag {wildcards.mag:q} --input {input.fasta:q} \
          --removals {input.removals:q} --trims {input.trims:q} --repairs {input.repairs:q} \
          --output {output.fasta:q} --report {output.report:q} > {log:q} 2>&1
        """


rule assess_deposit_mag_checkm2:
    input:
        fasta=f"{DEP}/cleaned_mags/{{mag}}.fa",
        report=f"{DEP}/edits/{{mag}}.json",
        database=DEPOSIT_CONFIG["databases"]["checkm2"],
    output:
        f"{DEP}/checkm2/{{mag}}/quality_report.tsv",
    params:
        tools=DEPOSIT_CONFIG["tools"]["checkm2"],
        outdir=lambda w: f"{DEP}/checkm2/{w.mag}",
        scratch=lambda w: f"tmp/deposit-c2/{list(DEPOSIT_MAGS).index(w.mag) + 1:03d}",
    wildcard_constraints:
        mag="|".join(sorted(DEPOSIT_EDITED)),
    threads: 8
    resources: mem_mb=32000, runtime=240, partition="ycga"
    log: f"{DEP}/logs/checkm2/{{mag}}.log"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        mkdir -p {params.scratch:q}
        export TMPDIR=$(realpath {params.scratch:q})
        checkm2 --version > {log:q} 2>&1
        checkm2 predict --threads {threads} --input {input.fasta:q} --extension fa \
          --database_path {input.database:q} --output-directory {params.outdir:q} --force >> {log:q} 2>&1
        """


rule package_deposit_sequences:
    input:
        code=["workflow/scripts/package_deposit_sequences.py", "workflow/scripts/prepare_deposit_mag.py"],
        config="config/deposit_preparation.json",
        decisions=DEPOSIT_CONFIG["decisions"],
        names=DEPOSIT_CONFIG["organism_names"],
        catalog=DEPOSIT_CONFIG["mag_catalog"],
        metadata=DEPOSIT_CONFIG["sample_manifest"],
        viral_classes=DEPOSIT_CONFIG["viral_classifications"],
        viral_endogenous_manifest=DEPOSIT_CONFIG["viral_endogenous_manifest"],
        fasta=expand(f"{DEP}/cleaned_mags/{{mag}}.fa", mag=sorted(DEPOSIT_MAGS)),
        edits=expand(f"{DEP}/edits/{{mag}}.json", mag=sorted(DEPOSIT_MAGS)),
        quality=expand(f"{DEP}/checkm2/{{mag}}/quality_report.tsv", mag=sorted(DEPOSIT_EDITED)),
        other=[path for spec in DEPOSIT_CONFIG["other_sequences"].values() for path in spec.values()],
    output:
        manifest=f"{DEP}/bundle/mag_manifest.tsv",
        quality=f"{DEP}/bundle/curated_mag_catalog.tsv",
        status=f"{DEP}/bundle/status.json",
        sets=expand(f"{DEP}/bundle/sets/{{seqset}}.fna", seqset=DEPOSIT_SET_IDS),
        mags=expand(f"{DEP}/bundle/mags/{{mag}}.fna", mag=sorted(DEPOSIT_MAGS)),
        metadata=f"{DEP}/bundle/mimag_biosample_metadata.tsv",
        mapping=f"{DEP}/bundle/sequence_id_mapping.tsv",
        set_manifest=f"{DEP}/bundle/sequence_sets.tsv",
        genbank=f"{DEP}/bundle/genbank_candidates.tsv",
        archive=f"{DEP}/bundle/archive_genomes.tsv",
        held=f"{DEP}/bundle/held_genomes.tsv",
    threads: 1
    resources: mem_mb=4000, runtime=30, partition="ycga"
    log: f"{DEP}/logs/package.log"
    shell:
        """
        {DEP_PY:q} workflow/scripts/package_deposit_sequences.py --config {input.config:q} > {log:q} 2>&1
        """


rule prepare_deposit_validation_inputs:
    input:
        code="workflow/scripts/prepare_deposit_validation_inputs.py",
        package_code="workflow/scripts/package_deposit_sequences.py",
        manifest=rules.package_deposit_sequences.output.manifest,
        metadata=rules.package_deposit_sequences.output.metadata,
        sets=rules.package_deposit_sequences.output.set_manifest,
        sequences=rules.package_deposit_sequences.output.sets,
        mags=rules.package_deposit_sequences.output.mags,
    output:
        directory(f"{DEP}/bundle/validation_inputs"),
    params: bundle=f"{DEP}/bundle"
    threads: 1
    resources: mem_mb=2000, runtime=15, partition="ycga"
    log: f"{DEP}/logs/validation_inputs.log"
    shell:
        "{DEP_PY:q} {input.code:q} --bundle {params.bundle:q} > {log:q} 2>&1"


rule screen_deposit_adaptors:
    input:
        fasta=f"{DEP}/bundle/sets/{{seqset}}.fna",
        runner=DEPOSIT_CONFIG["tools"]["fcs_adaptor"],
        image=DEPOSIT_CONFIG["tools"]["fcs_adaptor_image"],
    output:
        report=f"{DEP}/fcs_adaptor/{{seqset}}/fcs_adaptor_report.txt",
    params:
        outdir=lambda w: f"{DEP}/fcs_adaptor/{w.seqset}",
        tax=lambda w: "--prok" if w.seqset in ("mags", "viruses", "viral_endogenous") else "--euk",
        scratch=lambda w: f"{DEP}/tmp/fcs_adaptor/{w.seqset}",
    wildcard_constraints:
        seqset="|".join(DEPOSIT_SET_IDS),
    threads: 8
    resources: mem_mb=32000, runtime=240, partition="ycga"
    log: f"{DEP}/logs/fcs_adaptor/{{seqset}}.log"
    shell:
        """
        mkdir -p {params.scratch:q} {params.outdir:q}
        export TMPDIR=$(realpath {params.scratch:q})
        export APPTAINER_TMPDIR="$TMPDIR" SINGULARITY_TMPDIR="$TMPDIR"
        bash {input.runner:q} --fasta-input {input.fasta:q} --output-dir {params.outdir:q} \
          {params.tax} --container-engine singularity --image {input.image:q} > {log:q} 2>&1
        """


rule check_deposit_fcsgx_inputs:
    input:
        code="workflow/scripts/check_deposit_fcsgx_inputs.py",
        helpers=["workflow/scripts/package_deposit_sequences.py", "workflow/scripts/prepare_deposit_mag.py"],
        config="config/deposit_preparation.json",
        manifest=rules.package_deposit_sequences.output.manifest,
        validation_drafts=rules.prepare_deposit_validation_inputs.output,
        fasta=rules.package_deposit_sequences.output.mags,
        names=DEPOSIT_CONFIG["organism_names"],
        divisions=DEPOSIT_CONFIG["fcs_gx_division_overrides"],
        database=f"{DEPOSIT_CONFIG['databases']['fcs_gx']}/all.blast_div.tsv.gz",
        adaptor=f"{DEP}/fcs_adaptor/mags/fcs_adaptor_report.txt",
    output:
        f"{DEP}/fcs_gx_preflight.tsv",
    threads: 1
    resources: mem_mb=2000, runtime=15, partition="ycga"
    log: f"{DEP}/logs/fcs_gx_preflight.log"
    shell:
        "{DEP_PY:q} {input.code:q} --config {input.config:q} --output {output:q} > {log:q} 2>&1"


rule screen_deposit_mag_fcsgx:
    input:
        fasta=f"{DEP}/bundle/mags/{{mag}}.fna",
        preflight=rules.check_deposit_fcsgx_inputs.output,
        names=DEPOSIT_CONFIG["organism_names"],
        divisions=DEPOSIT_CONFIG["fcs_gx_division_overrides"],
        adaptor=f"{DEP}/fcs_adaptor/mags/fcs_adaptor_report.txt",
        runner=DEPOSIT_CONFIG["tools"]["fcs_gx"],
        image=DEPOSIT_CONFIG["tools"]["fcs_gx_image"],
    output:
        directory(f"{DEP}/fcs_gx/{{mag}}"),
    params:
        database=DEPOSIT_CONFIG["databases"]["fcs_gx"],
        taxid=deposit_taxid,
        division=deposit_division_options,
        enabled=lambda w: str(config.get("enable_fcsgx", False)).lower(),
        scratch=lambda w: f"{DEP}/tmp/fcs_gx/{w.mag}",
    threads: 16
    resources: mem_mb=650000, runtime=1440, partition="bigmem"
    log: f"{DEP}/logs/fcs_gx/{{mag}}.log"
    shell:
        """
        test {params.enabled:q} = true || {{ echo 'FCS-GX requires enable_fcsgx=true.' >&2; exit 2; }}
        mkdir -p {params.scratch:q} {output:q}
        export TMPDIR=$(realpath {params.scratch:q})
        export APPTAINER_TMPDIR="$TMPDIR" SINGULARITY_TMPDIR="$TMPDIR"
        export FCS_DEFAULT_IMAGE={input.image:q} GX_NUM_CORES={threads}
        {DEP_PY:q} {input.runner:q} screen genome --fasta {input.fasta:q} --out-dir {output:q} \
          --gx-db {params.database:q} --tax-id {params.taxid} {params.division:q} > {log:q} 2>&1
        """


rule deposit_preparation:
    input:
        rules.package_deposit_sequences.output,
        rules.prepare_deposit_validation_inputs.output,
        expand(f"{DEP}/fcs_adaptor/{{seqset}}/fcs_adaptor_report.txt", seqset=DEPOSIT_SET_IDS),


rule deposit_fcsgx:
    input:
        expand(f"{DEP}/fcs_gx/{{mag}}", mag=sorted(DEPOSIT_MAGS)),
