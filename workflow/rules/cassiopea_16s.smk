"""
Do the siphonophore Mycoplasmatales, or close relatives, occur in the upside-down
jellyfish Cassiopea? Florida Keys Cassiopea 16S amplicons (external and
gastrovascular tissue, with water, sediment and swab controls) are denoised and
matched to the Mycoplasmatales MAG 16S genes (config/cassiopea_16s.json; run
manifest with MD5s in config/cassiopea_16s_runs.tsv).

1. Obtain reads from ENA with MD5 checks: download_cassiopea_reads.
2. Merge, primer-trim and filter Illumina pairs, then denoise:
   merge_amplicon_pairs_vsearch, denoise_amplicons_vsearch.
3. Match denoised variants and PacBio full-length reads to the MAG 16S genes:
   match_cassiopea_16s_vsearch.
4. Summarize matches by sample type: summarize_cassiopea_matches.
"""

CAS = f"{WORK}/cassiopea_16s"
with open(ROOT / "config" / "cassiopea_16s.json") as handle:
    CAS_CONFIG = json.load(handle)
CAS_RUNS_PATH = str(ROOT / CAS_CONFIG["runs"])
CAS_TOOLS = f"{CAS_CONFIG['tool_prefixes']['vsearch']}/bin"


def cassiopea_resources(name):
    return {key: CAS_CONFIG["resources"][name][key] for key in ("mem_mb", "runtime", "partition")}


rule download_cassiopea_reads:
    input:
        runs=CAS_RUNS_PATH,
    output:
        directory(f"{CAS}/reads"),
    log:
        f"{WORK}/logs/cassiopea_16s/download.log",
    threads: 1
    resources:
        **cassiopea_resources("download"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        mkdir -p {output:q}
        tail -n +2 {input.runs:q} | cut -f 11,12 | while IFS=$'\\t' read -r urls md5s; do
          IFS=';' read -ra url <<< "$urls"
          IFS=';' read -ra md5 <<< "$md5s"
          for i in "${{!url[@]}}"; do
            file={output:q}/$(basename "${{url[$i]}}")
            curl -fsSL --retry 5 -o "$file" "${{url[$i]}}"
            echo "${{md5[$i]}}  $file" | md5sum -c --quiet || {{ echo "MD5 mismatch: $file" >&2; exit 1; }}
          done
        done > {log:q} 2>&1
        """


rule merge_amplicon_pairs_vsearch:
    input:
        runs=CAS_RUNS_PATH,
        reads=f"{CAS}/reads",
    output:
        f"{CAS}/illumina_filtered.fasta",
    params:
        tools=CAS_TOOLS,
        max_diffs=CAS_CONFIG["illumina"]["merge_max_diffs"],
        max_ee=CAS_CONFIG["illumina"]["max_expected_errors"],
        min_length=CAS_CONFIG["illumina"]["min_length"],
        max_length=CAS_CONFIG["illumina"]["max_length"],
        strip_left=len(CAS_CONFIG["primers"]["forward"]),
        strip_right=len(CAS_CONFIG["primers"]["reverse"]),
    log:
        f"{WORK}/logs/cassiopea_16s/merge.log",
    threads: CAS_CONFIG["resources"]["amplicons"]["threads"]
    resources:
        **cassiopea_resources("amplicons"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        vsearch --version > {log:q} 2>&1
        tmp=$(mktemp -d {output:q}.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        # Labels start with the run accession, which otutab reads as the sample.
        awk -F '\\t' 'NR > 1 && $4 == "ILLUMINA" {{print $1}}' {input.runs:q} | while read -r run; do
          vsearch --fastq_mergepairs {input.reads:q}/"$run"_1.fastq.gz --reverse {input.reads:q}/"$run"_2.fastq.gz \
            --fastq_maxdiffs {params.max_diffs} --threads {threads} --fastqout "$tmp/merged.fq" >> {log:q} 2>&1
          vsearch --fastx_filter "$tmp/merged.fq" --fastq_stripleft {params.strip_left} \
            --fastq_stripright {params.strip_right} --fastq_maxee {params.max_ee} \
            --fastq_minlen {params.min_length} --fastq_maxlen {params.max_length} --fastq_maxns 0 \
            --relabel "$run." --fastaout "$tmp/$run.fasta" >> {log:q} 2>&1
        done
        cat "$tmp"/SRR*.fasta > {output:q}.tmp
        mv {output:q}.tmp {output:q}
        """


rule denoise_amplicons_vsearch:
    input:
        reads=f"{CAS}/illumina_filtered.fasta",
    output:
        zotus=f"{CAS}/zotus.fasta",
        table=f"{CAS}/zotu_table.tsv",
    params:
        tools=CAS_TOOLS,
        min_size=CAS_CONFIG["illumina"]["unoise_min_size"],
        identity=CAS_CONFIG["illumina"]["otutab_identity"],
    log:
        f"{WORK}/logs/cassiopea_16s/denoise.log",
    threads: CAS_CONFIG["resources"]["amplicons"]["threads"]
    resources:
        **cassiopea_resources("amplicons"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        tmp=$(mktemp -d {output.zotus:q}.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        vsearch --derep_fulllength {input.reads:q} --sizeout --minuniquesize 2 --output "$tmp/uniques.fa" > {log:q} 2>&1
        vsearch --cluster_unoise "$tmp/uniques.fa" --minsize {params.min_size} --centroids "$tmp/denoised.fa" \
          --threads {threads} >> {log:q} 2>&1
        vsearch --uchime3_denovo "$tmp/denoised.fa" --nonchimeras {output.zotus:q}.tmp --relabel ZOTU \
          >> {log:q} 2>&1
        vsearch --usearch_global {input.reads:q} --db {output.zotus:q}.tmp --id {params.identity} --strand plus \
          --threads {threads} --otutabout {output.table:q}.tmp >> {log:q} 2>&1
        mv {output.zotus:q}.tmp {output.zotus:q}
        mv {output.table:q}.tmp {output.table:q}
        """


rule match_cassiopea_16s_vsearch:
    input:
        runs=CAS_RUNS_PATH,
        reads=f"{CAS}/reads",
        zotus=f"{CAS}/zotus.fasta",
        genes=f"{WORK}/cnidarian_16s/mycoplasmatales_mag_16s.fasta",
    output:
        zotus=f"{CAS}/zotu_matches.tsv",
        pacbio=f"{CAS}/pacbio_matches.tsv",
        pacbio_reads=f"{CAS}/pacbio_reads.tsv",
    params:
        tools=CAS_TOOLS,
        min_identity=CAS_CONFIG["match"]["min_identity"],
        coverage=CAS_CONFIG["match"]["query_coverage"],
        min_length=CAS_CONFIG["pacbio"]["min_length"],
        max_length=CAS_CONFIG["pacbio"]["max_length"],
    log:
        f"{WORK}/logs/cassiopea_16s/match.log",
    threads: CAS_CONFIG["resources"]["amplicons"]["threads"]
    resources:
        **cassiopea_resources("amplicons"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        export PATH={params.tools:q}:$PATH
        tmp=$(mktemp -d {output.zotus:q}.XXXXXX)
        trap 'rm -rf "$tmp"' EXIT
        # Denoised variants must align over their full length; best MAG gene per variant.
        vsearch --usearch_global {input.zotus:q} --db {input.genes:q} --id {params.min_identity} \
          --query_cov {params.coverage} --strand both --maxaccepts 0 --maxrejects 0 --threads {threads} \
          --userfields query+target+id+qcov --userout "$tmp/zotus.tsv" > {log:q} 2>&1
        printf 'zotu\\tmag_16s\\tidentity\\tcoverage\\n' > {output.zotus:q}.tmp
        sort -t $'\\t' -k1,1 -k3,3gr "$tmp/zotus.tsv" | awk -F '\\t' -v OFS='\\t' '!seen[$1]++ {{print $1, $2, $3/100, $4/100}}' \
          >> {output.zotus:q}.tmp
        # PacBio HiFi full-length reads: count per run, then best MAG gene per read.
        printf 'run\\treads_in_length_range\\n' > {output.pacbio_reads:q}.tmp
        printf 'run\\tread\\tmag_16s\\tidentity\\tcoverage\\n' > {output.pacbio:q}.tmp
        awk -F '\\t' 'NR > 1 && $4 == "PACBIO_SMRT" {{print $1}}' {input.runs:q} | while read -r run; do
          # HiFi reads carry quality values up to 93.
          vsearch --fastx_filter {input.reads:q}/"$run"_subreads.fastq.gz --fastq_qmax 93 --fastq_minlen {params.min_length} \
            --fastq_maxlen {params.max_length} --fastaout "$tmp/$run.fa" >> {log:q} 2>&1
          printf '%s\\t%s\\n' "$run" "$(grep -c '>' "$tmp/$run.fa")" >> {output.pacbio_reads:q}.tmp
          vsearch --usearch_global "$tmp/$run.fa" --db {input.genes:q} --id {params.min_identity} \
            --strand both --maxaccepts 0 --maxrejects 0 --threads {threads} \
            --userfields query+target+id+qcov --userout "$tmp/$run.tsv" >> {log:q} 2>&1
          sort -t $'\\t' -k1,1 -k3,3gr "$tmp/$run.tsv" \
            | awk -F '\\t' -v OFS='\\t' -v run="$run" '!seen[$1]++ {{print run, $1, $2, $3/100, $4/100}}' \
            >> {output.pacbio:q}.tmp
        done
        mv {output.zotus:q}.tmp {output.zotus:q}
        mv {output.pacbio:q}.tmp {output.pacbio:q}
        mv {output.pacbio_reads:q}.tmp {output.pacbio_reads:q}
        """


rule summarize_cassiopea_matches:
    input:
        script="workflow/scripts/summarize_cassiopea_matches.py",
        runs=CAS_RUNS_PATH,
        table=f"{CAS}/zotu_table.tsv",
        zotus=f"{CAS}/zotu_matches.tsv",
        pacbio=f"{CAS}/pacbio_matches.tsv",
        pacbio_reads=f"{CAS}/pacbio_reads.tsv",
    output:
        samples=f"{CAS}/summary/sample_matches.tsv",
        types=f"{CAS}/summary/sample_type_summary.tsv",
        pacbio=f"{CAS}/summary/pacbio_summary.tsv",
    params:
        config=str(ROOT / "config" / "cassiopea_16s.json"),
    log:
        f"{WORK}/logs/cassiopea_16s/summarize.log",
    threads: 1
    resources:
        **cassiopea_resources("summary"),
    conda:
        "../../envs/workflow.yaml"
    shell:
        """
        python {input.script:q} --config {params.config:q} --runs {input.runs:q} --table {input.table:q} \
          --zotus {input.zotus:q} --pacbio {input.pacbio:q} --pacbio-reads {input.pacbio_reads:q} \
          --samples {output.samples:q} --types {output.types:q} --pacbio-summary {output.pacbio:q} > {log:q} 2>&1
        """


rule cassiopea_16s:
    input:
        f"{CAS}/summary/sample_matches.tsv",
        f"{CAS}/summary/sample_type_summary.tsv",
        f"{CAS}/summary/pacbio_summary.tsv",
