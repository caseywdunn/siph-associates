#!/usr/bin/env Rscript
# Physalia analyses at a selected breadth threshold. Physical flowcells and
# inferential eligibility come from phase6_metadata.py, shared with batch tests.
suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(lme4); library(vegan)
})
args <- commandArgs(trailingOnly = TRUE)
analysis <- fromJSON(args[1])$physalia_models
all_grades <- fread(args[2]); libraries <- fread(args[3]); metadata <- fread(args[4])
out_design <- args[5]; out_models <- args[6]; out_permanova <- args[7]; out_permutation <- args[8]
grade_column <- if (length(args) >= 9) args[9] else "grade"
cores <- if (length(args) >= 10) as.integer(args[10]) else 1L
stopifnot(grade_column %in% c("grade", "grade_at_5pct", "grade_at_20pct"), cores >= 1)
set.seed(20260929)

# Multi-flowcell libraries remain in descriptive summaries but not inference.
phys <- merge(metadata[grepl("^Physalia", host_species) &
                         tolower(as.character(inference_eligible)) == "true" & ocean_region != "",
                       .(sample_id, ocean_region, flowcell)],
              libraries[, .(sample_id, input_pairs)], by = "sample_id")
stopifnot(!anyDuplicated(phys$sample_id), all(phys$input_pairs > 0), all(phys$flowcell != ""))
phys[, log_pairs := log10(input_pairs)]
phys[, ocean_region := factor(ocean_region)]

identifiable_design <- function(d) {
  design <- d[, .N, by = .(ocean_region, flowcell)]
  regions_ok <- design[, .(flowcells = uniqueN(flowcell)), by = ocean_region][flowcells >= 2, .N] >= 2
  spanning <- design[, .(regions = uniqueN(ocean_region)), by = flowcell][regions > 1, .N] > 0
  regions_ok && spanning
}
design <- phys[, .N, by = .(ocean_region, flowcell)]
design[, `:=`(grade_column = grade_column, inference_policy = "single_flowcell")]
fwrite(design[order(ocean_region, flowcell)], out_design, sep = "\t")
identifiable <- identifiable_design(phys)
cat("Physalia eligible libraries with region:", nrow(phys), "identifiable:", identifiable,
    "grade column:", grade_column, "\n")

present <- all_grades[role != "decoy" & get(grade_column) %in% c("validated", "high_confidence") &
                        sample_id %in% phys$sample_id]
counts <- present[, .N, by = target_id][N >= analysis$min_presences][order(target_id)]
results <- list()
if (identifiable) {
  for (target in counts$target_id) {
    d <- copy(phys)
    d[, presence := as.integer(sample_id %in% present[target_id == target, sample_id])]
    fit <- function(f) tryCatch(glmer(f, data = d, family = binomial,
                                      control = glmerControl(optimizer = "bobyqa")),
                                error = function(e) NULL)
    full <- fit(presence ~ log_pairs + ocean_region + (1 | flowcell))
    null <- fit(presence ~ log_pairs + (1 | flowcell))
    if (is.null(full) || is.null(null)) {
      results[[target]] <- data.table(target_id = target, presences = sum(d$presence), chisq = NA_real_,
                                      df = NA_real_, p_value = NA_real_, status = "fit_failed", messages = "")
    } else {
      test <- anova(null, full)
      messages <- c(full@optinfo$conv$lme4$messages, null@optinfo$conv$lme4$messages)
      status <- if (length(messages)) "convergence_warning" else if (isSingular(full) || isSingular(null))
        "fitted_singular" else "fitted"
      results[[target]] <- data.table(target_id = target, presences = sum(d$presence),
                                      chisq = test$Chisq[2], df = test$Df[2], p_value = test$`Pr(>Chisq)`[2],
                                      status = status,
                                      messages = gsub("[\t\n]+", " ", paste(unique(messages), collapse = "; ")))
    }
  }
}
models <- if (length(results)) rbindlist(results) else
  data.table(target_id = character(), presences = integer(), chisq = numeric(), df = numeric(),
             p_value = numeric(), status = character(), messages = character())
if (!identifiable) models <- data.table(target_id = counts$target_id, presences = counts$N,
                                       chisq = NA_real_, df = NA_real_, p_value = NA_real_,
                                       status = "region_not_identifiable", messages = "")
models[, `:=`(q_value = p.adjust(p_value, method = "BH"), libraries = nrow(phys),
               grade_column = grade_column)]
fwrite(models, out_models, sep = "\t")

# Approved replacement statistic: logistic deviance gain beyond read depth.
# Cache model matrices and the null fit rather than parsing 10,000 formulas per taxon.
# The same seeded within-flowcell label permutations are used for every target.
permutation <- data.table(target_id = character(), presences = integer(), deviance_gain = numeric(),
                          permutations = integer(), p_value = numeric(), status = character())
if (identifiable && nrow(counts)) {
  null_x <- model.matrix(~ log_pairs, data = phys)
  region_x <- model.matrix(~ ocean_region, data = phys)[, -1, drop = FALSE]
  set.seed(analysis$deviation_2026_09_30$seed)
  blocks <- split(seq_len(nrow(phys)), phys$flowcell)
  permutations <- replicate(analysis$permutations, {
    indices <- seq_len(nrow(phys))
    for (idx in blocks) if (length(idx) > 1) indices[idx] <- idx[sample.int(length(idx))]
    indices
  })
  fit_deviance <- function(x, y) suppressWarnings(glm.fit(x, y, family = binomial())$deviance)
  test_target <- function(target) {
    y <- as.integer(phys$sample_id %in% present[target_id == target, sample_id])
    null_deviance <- fit_deviance(null_x, y)
    observed <- null_deviance - fit_deviance(cbind(null_x, region_x), y)
    null_gains <- vapply(seq_len(ncol(permutations)), function(i)
      null_deviance - fit_deviance(cbind(null_x, region_x[permutations[, i], , drop = FALSE]), y), numeric(1))
    if (!is.finite(observed) || any(!is.finite(null_gains))) stop("Non-finite permutation statistic: ", target)
    data.table(target_id = target, presences = sum(y), deviance_gain = observed,
               permutations = analysis$permutations,
               p_value = (sum(null_gains >= observed - 1e-9) + 1) / (analysis$permutations + 1),
               status = "tested")
  }
  permutation <- rbindlist(parallel::mclapply(counts$target_id, test_target,
                                             mc.cores = min(cores, nrow(counts)), mc.set.seed = FALSE))
}
if (!identifiable && nrow(counts)) permutation <- data.table(
  target_id = counts$target_id, presences = counts$N, deviance_gain = NA_real_,
  permutations = 0L, p_value = NA_real_, status = "region_not_identifiable")
permutation[, `:=`(q_value = p.adjust(p_value, method = "BH"), libraries = nrow(phys),
                    grade_column = grade_column)]
fwrite(permutation, out_permutation, sep = "\t")

# Jaccard distances among non-empty communities only. Keep the excluded count
# explicit: a double-empty pair has no defined binary Jaccard denominator.
matrix_wide <- dcast(present[, .(sample_id, target_id, value = 1L)], sample_id ~ target_id,
                     value.var = "value", fill = 0L)
d <- merge(phys, matrix_wide, by = "sample_id", all.x = TRUE)
target_columns <- setdiff(names(matrix_wide), "sample_id")
community <- as.matrix(d[, ..target_columns])
community[is.na(community)] <- 0L
keep <- rowSums(community) > 0
community_identifiable <- sum(keep) > 2 && identifiable_design(d[keep])
if (community_identifiable) {
  set.seed(20260929)
  perm <- how(nperm = analysis$permutations, blocks = factor(d$flowcell[keep]))
  fit <- adonis2(vegdist(community[keep, , drop = FALSE], method = "jaccard", binary = TRUE) ~
                   log_pairs + ocean_region, data = droplevels(d[keep]), permutations = perm, by = "terms")
  permanova <- data.table(term = rownames(fit), as.data.table(fit))
  permanova[term == "log_pairs", term := "log10(input_pairs)"]
  permanova[, status := "tested"]
} else {
  permanova <- data.table(term = "ocean_region", Df = NA_real_, SumOfSqs = NA_real_, R2 = NA_real_,
                          F = NA_real_, `Pr(>F)` = NA_real_, status = "region_not_identifiable")
}
permanova[, `:=`(libraries = sum(keep), eligible_libraries = nrow(phys),
                 excluded_empty_communities = sum(!keep), region_identifiable = community_identifiable,
                 grade_column = grade_column)]
fwrite(permanova, out_permanova, sep = "\t")
cat("models:", nrow(models), "permanova libraries:", sum(keep), "\n")
