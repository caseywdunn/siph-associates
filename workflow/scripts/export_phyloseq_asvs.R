#!/usr/bin/env Rscript
# Export a McCauley et al. (2023) cnidarian 16S library (phyloseq RDS) as ASV FASTA and taxonomy.
# Usage: export_phyloseq_asvs.R <library.rdata> <library label> <out.fasta> <out.taxonomy.tsv>
suppressPackageStartupMessages({library(phyloseq); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
x <- readRDS(args[1]); label <- args[2]
ids <- paste0(label, "|", taxa_names(x))
seqs <- as.character(refseq(x))
writeLines(paste0(">", ids, "\n", seqs), args[3])
tax <- as.data.table(as(tax_table(x), "matrix"), keep.rownames = "asv")
tax[, asv := paste0(label, "|", asv)]
fwrite(tax, args[4], sep = "\t")
cat(label, ": ", ntaxa(x), " ASVs, ", nsamples(x), " samples\n", sep = "")
