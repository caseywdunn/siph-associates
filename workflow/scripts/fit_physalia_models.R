#!/usr/bin/env Rscript
# Primary analysis B2: Physalia presence models and PERMANOVA, fitted only when
# ocean region is identifiable apart from flowcell (docs/phase6_analysis_plan.md).
suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(lme4); library(vegan)
})
args <- commandArgs(trailingOnly = TRUE)
analysis <- fromJSON(args[1])$physalia_models
grades <- fread(args[2]); libraries <- fread(args[3]); manifest <- fread(args[4])
out_design <- args[5]; out_models <- args[6]; out_permanova <- args[7]
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
                                      df = NA, p_value = NA, status = "fit_failed")
    } else {
      test <- anova(null, full)
      singular <- isSingular(full)
      results[[target]] <- data.table(target_id = target, presences = sum(d$presence),
                                      chisq = test$Chisq[2], df = test$Df[2], p_value = test$`Pr(>Chisq)`[2],
                                      status = if (singular) "fitted_singular" else "fitted")
    }
  }
}
models <- if (length(results)) rbindlist(results) else
  data.table(target_id = character(), presences = integer(), chisq = numeric(), df = integer(),
             p_value = numeric(), status = character())
if (!identifiable) models <- data.table(target_id = counts$target_id, presences = counts$N, chisq = NA,
                                        df = NA, p_value = NA, status = "region_not_identifiable")
models[, q_value := p.adjust(p_value, method = "BH")]
fwrite(models, out_models, sep = "\t")

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
