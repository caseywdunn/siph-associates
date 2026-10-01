"""
Where across Cnidaria do 16S variants matching the siphonophore Mycoplasmatales
occur? Full-length 16S genes from the Mycoplasmatales MAGs are compared with the
unified Cnidarian Microbiome Database of McCauley et al. (2023)
(config/cnidarian_16s.json).

1. Obtain the database libraries with MD5 checks: fetch_cnidarian_16s_library.
2. Export each library's ASVs: export_cnidarian_asvs_phyloseq.
3. Extract 16S genes from the Mycoplasmatales MAGs: extract_mycoplasmatales_16s.
4. Match ASVs over their full length to the MAG 16S genes: match_asvs_mycoplasmatales_vsearch.
5. Summarize the hosts carrying matching ASVs: summarize_matching_asv_hosts.
"""

C16 = f"{WORK}/cnidarian_16s"
with open(ROOT / "config" / "cnidarian_16s.json") as handle:
    C16_CONFIG = json.load(handle)
C16_LIBRARIES = sorted(C16_CONFIG["database"]["files"])

wildcard_constraints:
    library="|".join(C16_LIBRARIES),


rule fetch_cnidarian_16s_library:
    output:
        f"{C16}/mccauley2023/{{library}}_Library.rdata",
    params:
        file_id=lambda wildcards: C16_CONFIG["database"]["files"][wildcards.library]["figshare_file"],
        md5=lambda wildcards: C16_CONFIG["database"]["files"][wildcards.library]["md5"],
    log:
        f"{WORK}/logs/cnidarian_16s/fetch.{{library}}.log",
    threads: 1
    resources:
        **phase4_resources("download"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        curl -fsSL --retry 5 -o {output:q}.part https://ndownloader.figshare.com/files/{params.file_id} 2> {log:q}
        echo "{params.md5}  {output}.part" | md5sum -c - >> {log:q} 2>&1
        mv {output:q}.part {output:q}
        """


rule export_cnidarian_asvs_phyloseq:
    input:
        script="workflow/scripts/export_phyloseq_asvs.R",
        library=f"{C16}/mccauley2023/{{library}}_Library.rdata",
    output:
        fasta=f"{C16}/asvs/{{library}}.fasta",
        taxonomy=f"{C16}/asvs/{{library}}.taxonomy.tsv",
    params:
        tools=f"{C16_CONFIG['tool_prefixes']['phyloseq']}/bin",
    log:
        f"{WORK}/logs/cnidarian_16s/export.{{library}}.log",
    threads: 1
    resources:
        **phase4_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        {params.tools}/Rscript {input.script:q} {input.library:q} {wildcards.library} \
          {output.fasta:q} {output.taxonomy:q} > {log:q} 2>&1
        """


rule extract_mycoplasmatales_16s:
    input:
        script="workflow/scripts/extract_mag_16s.py",
        mag_catalog=f"{WORK}/phase3_catalog/mags/mag_catalog.tsv",
    output:
        f"{C16}/mycoplasmatales_mag_16s.fasta",
    params:
        results=f"{WORK}/phase3_cohort",
        min_length=C16_CONFIG["mag_16s_min_length"],
    log:
        f"{WORK}/logs/cnidarian_16s/extract_mag_16s.log",
    threads: 1
    resources:
        **phase4_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --mag-catalog {input.mag_catalog:q} --results {params.results:q} \
          --order Mycoplasmatales --min-length {params.min_length} --output {output:q} > {log:q} 2>&1
        """


rule match_asvs_mycoplasmatales_vsearch:
    input:
        asvs=expand(f"{C16}/asvs/{{library}}.fasta", library=C16_LIBRARIES),
        genes=f"{C16}/mycoplasmatales_mag_16s.fasta",
    output:
        f"{C16}/asv_matches.tsv",
    params:
        tools=f"{C16_CONFIG['tool_prefixes']['vsearch']}/bin",
        min_identity=C16_CONFIG["match"]["min_identity"],
        coverage=C16_CONFIG["match"]["query_coverage"],
    log:
        f"{WORK}/logs/cnidarian_16s/match.log",
    threads: 8
    resources:
        **phase4_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        vsearch --version > {log:q} 2>&1
        printf 'library\\tasv\\tspecies_cluster\\tmag_16s\\tidentity\\tasv_coverage\\n' > {output:q}.tmp
        for asvs in {input.asvs:q}; do
          # ASVs are short; each must align over its full length to a MAG 16S gene.
          vsearch --usearch_global "$asvs" --db {input.genes:q} --id {params.min_identity} \
            --query_cov {params.coverage} --strand both --maxaccepts 0 --maxrejects 0 --threads {threads} \
            --userfields query+target+id+qcov --userout "$asvs.matches" >> {log:q} 2>&1
          awk -F '\\t' -v OFS='\\t' '{{split($1, a, "|"); split($2, b, "|"); print a[1], $1, b[1], $2, $3/100, $4/100}}' \
            "$asvs.matches" >> {output:q}.tmp
          rm "$asvs.matches"
        done
        mv {output:q}.tmp {output:q}
        """


rule summarize_matching_asv_hosts:
    input:
        script="workflow/scripts/summarize_asv_hosts.R",
        matches=f"{C16}/asv_matches.tsv",
        libraries=expand(f"{C16}/mccauley2023/{{library}}_Library.rdata", library=C16_LIBRARIES),
    output:
        samples=f"{C16}/matching_asv_samples.tsv",
        hosts=f"{C16}/matching_asv_hosts.tsv",
    params:
        tools=f"{C16_CONFIG['tool_prefixes']['phyloseq']}/bin",
    log:
        f"{WORK}/logs/cnidarian_16s/summarize.log",
    threads: 1
    resources:
        # Loading all six phyloseq libraries needs more than the 16 GB summary default.
        **{**phase4_resources("summary"), "mem_mb": 48000},
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        {params.tools}/Rscript {input.script:q} {input.matches:q} {output.samples:q} {output.hosts:q} \
          {input.libraries:q} > {log:q} 2>&1
        """


rule cnidarian_16s:
    input:
        f"{C16}/matching_asv_hosts.tsv",
