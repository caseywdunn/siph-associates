#!/usr/bin/env Rscript
# Host distribution of cnidarian-database ASVs that match siphonophore Mycoplasmatales 16S genes.
# Usage: summarize_asv_hosts.R <matches.tsv> <out.samples.tsv> <out.hosts.tsv> <library.rdata>...
suppressPackageStartupMessages({library(phyloseq); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
matches <- fread(args[1])
libraries <- args[-(1:3)]
rows <- list()
for (path in libraries) {
  label <- sub("_Library.rdata$", "", basename(path))
  hits <- matches[library == label]
  if (!nrow(hits)) next
  x <- readRDS(path)
  keep <- intersect(sub(".*\\|", "", hits$asv), taxa_names(x))
  if (!length(keep)) next
  counts <- as(otu_table(prune_taxa(keep, x)), "matrix")
  if (!taxa_are_rows(x)) counts <- t(counts)
  meta <- as.data.table(as(sample_data(x), "data.frame"), keep.rownames = "sample")
  long <- as.data.table(as.table(counts))
  setnames(long, c("asv", "sample", "reads"))
  long <- long[reads > 0]
  long[, asv := paste0(label, "|", asv)]
  long <- merge(long, hits[, .(asv, species_cluster, identity)], by = "asv")
  totals <- data.table(sample = colnames(counts), library_reads = sample_sums(x)[colnames(counts)])
  long <- merge(merge(long, totals, by = "sample"), meta, by = "sample")
  long[, `:=`(library = label, relative_abundance = reads / library_reads)]
  rows[[label]] <- long
}
samples <- rbindlist(rows, fill = TRUE)
fwrite(samples, args[2], sep = "\t")
hosts <- samples[, .(samples = uniqueN(sample), studies = uniqueN(Code),
                     median_relative_abundance = median(relative_abundance),
                     max_relative_abundance = max(relative_abundance),
                     best_identity = max(identity)),
                 by = .(species_cluster, Phylum, Class, Order, Family, Genus, Species)]
fwrite(hosts[order(species_cluster, -samples)], args[3], sep = "\t")
cat("matched samples:", uniqueN(samples$sample), " host species:", uniqueN(samples$Species), "\n")
