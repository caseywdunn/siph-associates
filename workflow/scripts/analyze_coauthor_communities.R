#!/usr/bin/env Rscript
# Conditional community comparisons of completed detections; no sequence access.
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(vegan)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
plots_only <- "--plots-only" %in% args[-1]
Sys.setenv(TZ = "UTC")
settings <- fromJSON(args[1], simplifyVector = FALSE)
out <- settings$output_dir
metadata <- fread(file.path(out, "specimen_context.tsv"))
if (!plots_only) {
membership <- fread(file.path(out, "contrast_membership.tsv"))
grades <- fread(settings$bacterial_grades)
if ("biological_include" %in% names(grades)) grades <- grades[tolower(as.character(biological_include)) == "true"]
stopifnot(!anyDuplicated(metadata$sample_id))
stopifnot(!anyDuplicated(grades[, .(sample_id, target_id)]))
stopifnot(all(grades$sample_id %in% metadata$sample_id))

# The number of distinct allocations is relevant when tissue groups are small;
# requesting 9,999 draws does not create more independent assignments.
label_allocations <- function(d) {
  blocks <- split(d$group, d$flowcell)
  exp(sum(vapply(blocks, function(labels) {
    counts <- table(labels)
    lgamma(sum(counts) + 1) - sum(lgamma(counts + 1))
  }, numeric(1))))
}

community_result <- function(spec, grade_column) {
  selected <- membership[contrast == spec$id & single_flowcell_eligible == TRUE]
  d <- merge(selected, metadata[, .(sample_id, log10_input_pairs)], by = "sample_id")
  d[, `:=`(group = factor(group), flowcell = factor(flowcell))]
  present <- grades[
    sample_id %in% d$sample_id & role != "decoy" &
      grepl("(^|;)d__Bacteria(;|$)", taxonomy) &
      get(grade_column) %in% c("validated", "high_confidence")
  ]
  taxa <- sort(unique(present$target_id))
  community <- matrix(0L, nrow(d), length(taxa),
    dimnames = list(d$sample_id, taxa)
  )
  if (nrow(present)) {
    community[cbind(match(present$sample_id, d$sample_id), match(present$target_id, taxa))] <- 1L
  }
  keep <- rowSums(community) > 0
  audit <- d[, .(contrast, sample_id, group, flowcell, input_pairs)]
  audit[, `:=`(grade_column = grade_column, bacterial_targets = rowSums(community),
               included_nonempty_profile = keep)]
  eligible <- nrow(d)
  d <- droplevels(d[keep])
  community <- community[keep, , drop = FALSE]
  possible <- label_allocations(d)
  result <- data.table(
    contrast = spec$id, grade_column = grade_column,
    eligible_libraries = eligible, included_libraries = nrow(d),
    excluded_empty_profiles = sum(!keep), groups = nlevels(d$group),
    physical_flowcells = nlevels(d$flowcell),
    distinct_label_allocations = possible,
    group_df = NA_real_, residual_df = NA_real_, r_squared = NA_real_,
    pseudo_f = NA_real_, permutations = 0L, p_value = NA_real_, status = "not_tested"
  )
  dispersion <- data.table(contrast = spec$id, grade_column = grade_column,
                          included_libraries = nrow(d), f_statistic = NA_real_,
                          p_value = NA_real_, status = "not_tested")
  finish <- function() list(test = result, dispersion = dispersion, membership = audit)
  if (nlevels(d$group) < 2 || any(table(d$group) < 3)) {
    result[, status := "insufficient_nonempty_profiles_per_group"]
    return(finish())
  }
  null_x <- model.matrix(~ log10_input_pairs + flowcell, data = d)
  full_x <- model.matrix(~ log10_input_pairs + flowcell + group, data = d)
  gained <- qr(full_x)$rank - qr(null_x)$rank
  residual_degrees <- nrow(d) - qr(full_x)$rank
  result[, `:=`(group_df = gained, residual_df = residual_degrees)]
  if (gained != nlevels(d$group) - 1) {
    result[, status := "group_confounded_with_flowcell_or_read_effort"]
    return(finish())
  }
  if (possible < 20 - 1e-6 || residual_degrees < 2) {
    result[, status := "insufficient_distinct_permutations_or_residual_df"]
    return(finish())
  }
  distances <- vegdist(community, method = "jaccard", binary = TRUE)
  permutation_design <- permute::how(nperm = settings$permutations, blocks = d$flowcell)
  set.seed(settings$seed)
  fit <- adonis2(distances ~ log10_input_pairs + flowcell + group,
    data = d, permutations = permutation_design, by = "margin"
  )
  stopifnot("group" %in% rownames(fit))
  result[, `:=`(
    r_squared = fit["group", "R2"], pseudo_f = fit["group", "F"],
    permutations = settings$permutations, p_value = fit["group", "Pr(>F)"], status = "tested"
  )]
  # Dispersion is outside the authorized revision; retain an explicit unrun status.
  dispersion[, status := "not_rerun_out_of_scope"]
  finish()
}

results <- list()
for (grade_column in unlist(settings$grade_columns)) {
  for (spec in settings$contrasts) {
    message(spec$id, ": ", grade_column)
    results[[paste(spec$id, grade_column)]] <- community_result(spec, grade_column)
  }
}
tests <- rbindlist(lapply(results, `[[`, "test"))
dispersion <- rbindlist(lapply(results, `[[`, "dispersion"))
members <- rbindlist(lapply(results, `[[`, "membership"))
# Correct the estimable contrasts within each threshold; sensitivities are
# reported separately from the primary 10% presence rule.
tests[, q_value := p.adjust(p_value, method = "BH"), by = grade_column]
fwrite(tests, file.path(out, "community_tests.tsv"), sep = "\t", na = "")
fwrite(dispersion, file.path(out, "dispersion_tests.tsv"), sep = "\t", na = "")
fwrite(members, file.path(out, "community_membership.tsv"), sep = "\t", na = "")
write_json(list(tests = tests, dispersion = dispersion,
                interpretation = settings$interpretation,
                model = settings$community_model),
  file.path(out, "community_values.json"), pretty = TRUE, auto_unbox = TRUE, na = "null"
)

}

theme_set(theme_classic(base_size = 10) + theme(
  strip.background = element_blank(), strip.text = element_text(face = "bold"),
  legend.position = "bottom", legend.box = "vertical",
  axis.text = element_text(colour = "black")
))

# Point and recorded interval values retain their original sampling uncertainty.
# The formalin specimen is identified by shape, without interpreting its high
# viral count as diversity comparable to the rest of the collection.
depth <- metadata[host_group == "non_Physalia"]
depth[, plot_group := fifelse(host_genus == "Nanomia", "Nanomia",
  fifelse(host_genus == "Bargmannia", "Bargmannia", "Other siphonophores"))]
depth[, preservation_flag := fifelse(background_control_flag, "Flagged formalin specimen", "Other specimens")]
depth_long <- melt(depth,
  id.vars = c("sample_id", "depth_m", "depth_min_m", "depth_max_m", "plot_group", "preservation_flag"),
  measure.vars = c("bacterial_targets", "viral_targets", "eukaryote_groups"),
  variable.name = "evidence", value.name = "detections"
)
depth_long[, evidence := factor(evidence,
  levels = c("bacterial_targets", "viral_targets", "eukaryote_groups"),
  labels = c("Bacterial genome targets", "Viral genome targets", "Eukaryote groups")
)]
depth_plot <- ggplot(depth_long, aes(x = depth_m, y = detections, colour = plot_group)) +
  geom_segment(aes(x = depth_min_m, xend = depth_max_m, yend = detections),
    alpha = 0.45, na.rm = TRUE
  ) +
  geom_point(aes(shape = preservation_flag), size = 2, alpha = 0.8) +
  facet_wrap(~ evidence, scales = "free_y", ncol = 1) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 10, 100, 300, 1000, 3000)) +
  scale_colour_manual(values = c("Bargmannia" = "#D55E00", "Nanomia" = "#0072B2", "Other siphonophores" = "#666666")) +
  scale_shape_manual(values = c("Flagged formalin specimen" = 17, "Other specimens" = 16)) +
  labs(x = "Collection depth (m)", y = "Accepted detections", colour = NULL, shape = NULL)
ggsave(file.path(out, "depth_patterns.pdf"), depth_plot, width = 7, height = 7.5, device = grDevices::cairo_pdf)
ggsave(file.path(out, "depth_patterns.png"), depth_plot, width = 7, height = 7.5, dpi = 180)

incidence <- fread(file.path(out, "targeted_incidence.tsv"))[
  scope == "all_libraries" & tissue == "tentacle"
]
incidence[, label := fifelse(subset == "utriculus_regions", ocean_region,
                            sub("Physalia ", "P. ", species_current))]
incidence[, panel := fifelse(subset == "utriculus_regions", "P. utriculus by region", "SW Pacific by species")]
incidence <- melt(incidence, id.vars = c("panel", "label", "total_libraries"),
  measure.vars = c("clade_a_positive_libraries", "clade_b_positive_libraries", "vibrio_positive_libraries"),
  variable.name = "associate", value.name = "positive"
)
incidence[, associate := factor(associate,
  levels = c("clade_a_positive_libraries", "clade_b_positive_libraries", "vibrio_positive_libraries"),
  labels = c("Clade A", "Clade B", "Vibrio")
)]
incidence_plot <- ggplot(incidence, aes(associate, label, fill = positive / total_libraries)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = paste0(positive, "/", total_libraries)), size = 3.2) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "white", high = "#80B1D3", limits = c(0, 1), name = "Fraction detected") +
  labs(x = NULL, y = NULL) + theme(axis.line = element_blank(), axis.ticks = element_blank())
ggsave(file.path(out, "targeted_incidence.pdf"), incidence_plot, width = 6, height = 6, device = grDevices::cairo_pdf)
ggsave(file.path(out, "targeted_incidence.png"), incidence_plot, width = 6, height = 6, dpi = 180)
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
if (!plots_only) print(tests)
