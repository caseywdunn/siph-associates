#!/usr/bin/env Rscript
# Primary analysis B2: Physalia presence models and PERMANOVA, fitted only when
# ocean region is identifiable apart from flowcell (docs/phase6_analysis_plan.md).
suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(lme4); library(vegan)
})
args <- commandArgs(trailingOnly = TRUE)
analysis <- fromJSON(args[1])$physalia_models
grades <- fread(args[2]); libraries <- fread(args[3]); manifest <- fread(args[4])
out_design <- args[5]; out_models <- args[6]; out_permanova <- args[7]; out_permutation <- args[8]
set.seed(20260929)

manifest <- manifest[tolower(include_primary) == "true"]
manifest[, sample_id := paste0(study, "__", sub(".*:", "", library_id))]
manifest[, flowcell := sapply(strsplit(sequencing_batches, ";"), function(b)
  paste(sort(unique(sub("^([^/]*/[^/]*/[^/]*).*", "\\1", b[b != ""]))), collapse = ";"))]
phys <- merge(manifest[grepl("^Physalia", species_current), .(sample_id, ocean_region, flowcell)],
              libraries[, .(sample_id, input_pairs)], by = "sample_id")
phys <- phys[ocean_region != ""]

# Identifiability: >= 2 regions each on >= 2 flowcells, and flowcells spanning regions.
design <- phys[, .N, by = .(ocean_region, flowcell)]
fwrite(design[order(ocean_region, flowcell)], out_design, sep = "\t")
regions_ok <- design[, .(flowcells = uniqueN(flowcell)), by = ocean_region][flowcells >= 2, .N] >= 2
spanning <- design[, .(regions = uniqueN(ocean_region)), by = flowcell][regions > 1, .N] > 0
identifiable <- regions_ok && spanning
cat("Physalia libraries with region:", nrow(phys), " identifiable:", identifiable, "\n")

present <- grades[role != "decoy" & grade %in% c("validated", "high_confidence") & sample_id %in% phys$sample_id]
counts <- present[, .N, by = target_id][N >= analysis$min_presences]
results <- list()
if (identifiable) {
  for (target in counts$target_id) {
    d <- copy(phys)
    d[, presence := as.integer(sample_id %in% present[target_id == target, sample_id])]
    d[, log_pairs := log10(input_pairs)]
    fit <- function(f) tryCatch(glmer(f, data = d, family = binomial,
                                      control = glmerControl(optimizer = "bobyqa")),
                                error = function(e) NULL)
    full <- fit(presence ~ log_pairs + ocean_region + (1 | flowcell))
    null <- fit(presence ~ log_pairs + (1 | flowcell))
    if (is.null(full) || is.null(null)) {
      results[[target]] <- data.table(target_id = target, presences = sum(d$presence), chisq = NA,
                                      df = NA, p_value = NA, status = "fit_failed", messages = "")
    } else {
      test <- anova(null, full)
      # Record singular fits and optimizer convergence messages so unreliable fits are visible.
      messages <- c(full@optinfo$conv$lme4$messages, null@optinfo$conv$lme4$messages)
      status <- if (length(messages)) "convergence_warning" else if (isSingular(full)) "fitted_singular" else "fitted"
      results[[target]] <- data.table(target_id = target, presences = sum(d$presence),
                                      chisq = test$Chisq[2], df = test$Df[2], p_value = test$`Pr(>Chisq)`[2],
                                      status = status,
                                      messages = gsub("[\t\n]+", " ", paste(unique(messages), collapse = "; ")))
    }
  }
}
models <- if (length(results)) rbindlist(results) else
  data.table(target_id = character(), presences = integer(), chisq = numeric(), df = integer(),
             p_value = numeric(), status = character(), messages = character())
if (!identifiable) models <- data.table(target_id = counts$target_id, presences = counts$N, chisq = NA,
                                        df = NA, p_value = NA, status = "region_not_identifiable",
                                        messages = "")
models[, q_value := p.adjust(p_value, method = "BH")]
fwrite(models, out_models, sep = "\t")

# Deviation (approved 2026-09-30): the pre-specified mixed models fail under separation, so the
# reported per-taxon region test is a permutation test mirroring the PERMANOVA design. The statistic
# is the deviance improvement from adding ocean region to a logistic model with log depth; region
# labels are permuted within flowcell, which stays valid when estimates do not converge.
deviance_gain <- function(presence, log_pairs, region) {
  null <- suppressWarnings(glm(presence ~ log_pairs, family = binomial))
  full <- suppressWarnings(glm(presence ~ log_pairs + region, family = binomial))
  deviance(null) - deviance(full)
}
permute_within <- function(labels, blocks) {
  out <- labels
  for (idx in split(seq_along(labels), blocks)) if (length(idx) > 1) out[idx] <- labels[idx][sample.int(length(idx))]
  out
}
permutation <- list()
if (identifiable) {
  d <- copy(phys)
  d[, log_pairs := log10(input_pairs)]
  for (target in counts$target_id) {
    set.seed(20260930)
    presence <- as.integer(d$sample_id %in% present[target_id == target, sample_id])
    observed <- deviance_gain(presence, d$log_pairs, factor(d$ocean_region))
    null_gains <- replicate(analysis$permutations,
                            deviance_gain(presence, d$log_pairs, factor(permute_within(d$ocean_region, d$flowcell))))
    permutation[[target]] <- data.table(target_id = target, presences = sum(presence),
                                        deviance_gain = observed, permutations = analysis$permutations,
                                        p_value = (sum(null_gains >= observed - 1e-9) + 1) / (analysis$permutations + 1))
  }
}
permutation <- if (length(permutation)) rbindlist(permutation) else
  data.table(target_id = character(), presences = integer(), deviance_gain = numeric(),
             permutations = integer(), p_value = numeric())
permutation[, q_value := p.adjust(p_value, method = "BH")]
fwrite(permutation, out_permutation, sep = "\t")

# Community PERMANOVA on validated presence, permutations restricted within flowcell.
matrix_wide <- dcast(present[, .(sample_id, target_id, value = 1L)], sample_id ~ target_id,
                     value.var = "value", fill = 0L)
d <- merge(phys, matrix_wide, by = "sample_id")
community <- as.matrix(d[, setdiff(names(matrix_wide), "sample_id"), with = FALSE])
keep <- rowSums(community) > 0
perm <- how(nperm = analysis$permutations, blocks = factor(d$flowcell[keep]))
fit <- adonis2(vegdist(community[keep, , drop = FALSE], method = "jaccard", binary = TRUE) ~
                 log10(input_pairs) + ocean_region, data = d[keep], permutations = perm, by = "terms")
permanova <- data.table(term = rownames(fit), as.data.table(fit))
permanova[, `:=`(libraries = sum(keep), region_identifiable = identifiable)]
fwrite(permanova, out_permanova, sep = "\t")
cat("models:", nrow(models), " permanova libraries:", sum(keep), "\n")
