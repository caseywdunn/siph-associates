#!/usr/bin/env Rscript
# A bounded omission sensitivity, first reproducing both original primary tests.
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(vegan)
})

args <- commandArgs(trailingOnly = TRUE)
Sys.setenv(TZ = "UTC")
settings <- fromJSON(args[1], simplifyVector = FALSE)
original <- fromJSON(settings$original_config, simplifyVector = FALSE)
metadata <- fread(settings$metadata)
membership <- fread(settings$membership)
grades <- fread(settings$bacterial_grades)
previous <- fread(settings$original_tests)
stopifnot(!anyDuplicated(metadata$sample_id))
stopifnot(!anyDuplicated(grades[, .(sample_id, target_id)]))
stopifnot(all(grades$sample_id %in% metadata$sample_id))
stopifnot(settings$bacterial_grades == original$bacterial_grades)
stopifnot(settings$grade_column == "grade")
stopifnot(all(unlist(settings$contrasts) %in%
  vapply(original$contrasts, `[[`, character(1), "id")))

community_result <- function(contrast_id, omit) {
  selected <- membership[contrast == contrast_id & single_flowcell_eligible == TRUE]
  d <- merge(selected, metadata[, .(sample_id, log10_input_pairs)], by = "sample_id")
  d[, `:=`(group = factor(group), flowcell = factor(flowcell))]
  present <- grades[
    sample_id %in% d$sample_id & role != "decoy" &
      grepl("(^|;)d__Bacteria(;|$)", taxonomy) &
      get(settings$grade_column) %in% c("validated", "high_confidence")
  ]
  if (omit) present <- present[target_id != settings$omit_target]
  taxa <- sort(unique(present$target_id))
  community <- matrix(0L, nrow(d), length(taxa),
    dimnames = list(d$sample_id, taxa)
  )
  if (nrow(present)) {
    community[cbind(match(present$sample_id, d$sample_id), match(present$target_id, taxa))] <- 1L
  }
  keep <- rowSums(community) > 0
  audit <- d[, .(contrast, sample_id, group, flowcell, input_pairs)]
  audit[, `:=`(omit_cutibacterium = omit, bacterial_targets = rowSums(community),
               included_nonempty_profile = keep)]
  eligible <- nrow(d)
  d <- droplevels(d[keep])
  community <- community[keep, , drop = FALSE]
  possible <- exp(sum(vapply(split(d$group, d$flowcell), function(labels) {
    counts <- table(labels)
    lgamma(sum(counts) + 1) - sum(lgamma(counts + 1))
  }, numeric(1))))
  null_x <- model.matrix(~ log10_input_pairs + flowcell, data = d)
  full_x <- model.matrix(~ log10_input_pairs + flowcell + group, data = d)
  gained <- qr(full_x)$rank - qr(null_x)$rank
  residual_degrees <- nrow(d) - qr(full_x)$rank
  # Do not silently replace a now-unidentifiable sensitivity with a partial test.
  stopifnot(nlevels(d$group) >= 2, all(table(d$group) >= 3))
  stopifnot(gained == nlevels(d$group) - 1, residual_degrees >= 2, possible >= 20)
  distances <- vegdist(community, method = "jaccard", binary = TRUE)
  permutation_design <- permute::how(nperm = original$permutations, blocks = d$flowcell)
  set.seed(original$seed)
  fit <- adonis2(distances ~ log10_input_pairs + flowcell + group,
    data = d, permutations = permutation_design, by = "margin"
  )
  set.seed(original$seed)
  dispersion <- permutest(betadisper(distances, d$group, type = "median", sqrt.dist = TRUE),
    permutations = permutation_design
  )
  result <- data.table(
    contrast = contrast_id, grade_column = settings$grade_column,
    omit_cutibacterium = omit, eligible_libraries = eligible,
    included_libraries = nrow(d), excluded_empty_profiles = sum(!keep),
    groups = nlevels(d$group), physical_flowcells = nlevels(d$flowcell),
    distinct_label_allocations = possible, group_df = gained,
    residual_df = residual_degrees, r_squared = fit["group", "R2"],
    pseudo_f = fit["group", "F"], permutations = original$permutations,
    p_value = fit["group", "Pr(>F)"], status = "tested",
    dispersion_f = dispersion$tab[1, "F"], dispersion_p = dispersion$tab[1, "Pr(>F)"]
  )
  if (!omit) {
    baseline <- previous[contrast == contrast_id & grade_column == settings$grade_column]
    stopifnot(nrow(baseline) == 1)
    for (field in c("eligible_libraries", "included_libraries", "excluded_empty_profiles",
                    "groups", "physical_flowcells", "group_df", "residual_df",
                    "r_squared", "pseudo_f", "permutations", "p_value")) {
      if (!isTRUE(all.equal(result[[field]], baseline[[field]], tolerance = 1e-10))) {
        stop("Baseline reconstruction differs: ", contrast_id, " / ", field)
      }
    }
  }
  list(test = result, membership = audit)
}

results <- list()
for (omit in c(FALSE, TRUE)) {
  for (contrast_id in unlist(settings$contrasts)) {
    message(contrast_id, ": omit BAC00019 = ", omit)
    results[[paste(contrast_id, omit)]] <- community_result(contrast_id, omit)
  }
}
tests <- rbindlist(lapply(results, `[[`, "test"))
tests[, q_value := p.adjust(p_value, method = "BH"), by = omit_cutibacterium]
members <- rbindlist(lapply(results, `[[`, "membership"))
dir.create(settings$output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(tests, file.path(settings$output_dir, "community_tests.tsv"), sep = "\t", na = "")
fwrite(members, file.path(settings$output_dir, "community_membership.tsv"), sep = "\t", na = "")
write_json(list(tests = tests, baseline_reconstruction_matches = TRUE,
                interpretation = settings$interpretation, model = original$community_model,
                omitted_target = settings$omit_target, permutations = original$permutations,
                seed = original$seed),
  file.path(settings$output_dir, "summary.json"), pretty = TRUE, auto_unbox = TRUE, na = "null"
)
values <- c("% Sensitivity excluding C. acnes BAC00019; baseline analysis retained.")
for (contrast_id in unlist(settings$contrasts)) {
  prefix <- if (contrast_id == "utriculus_pacific_regions") "cutibSensitivityWithinHost" else "cutibSensitivityWithinRegion"
  result <- tests[contrast == contrast_id & omit_cutibacterium == TRUE]
  suffixes <- c("Libraries", "RtwoPercent", "P", "Q")
  numbers <- c(result$included_libraries, sprintf("%.1f", 100 * result$r_squared),
               sprintf("%.4f", result$p_value), sprintf("%.4f", result$q_value))
  values <- c(values, paste0("\\newcommand{\\", prefix, suffixes, "}{", numbers, "}"))
}
writeLines(values, file.path(settings$output_dir, "values.tex"))
writeLines(capture.output(sessionInfo()), file.path(settings$output_dir, "session_info.txt"))
print(tests)
