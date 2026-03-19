################################################################################
## I-WABA Real Data Analysis: 14-Cancer Gene Expression Dataset
## Revised top-50 + all-gene screening workflow
##
## Outputs:
##   - cancer_top50_level1_v2.csv
##   - cancer_top50_level2_v2.csv
##   - cancer_top50_level2_aggregate_v2.csv
##   - cancer_top50_level3_v2.csv
##   - cancer_top50_level4_v2.csv
##   - cancer_top50_level5_v2.csv
##   - cancer_top50_summary_v2.csv
##   - cancer_top50_screening_v1.csv
##   - cancer_top50_inference_v1.csv
##   - cancer_allgenes_summary_v2.csv
##   - cancer_allgenes_screening_v1.csv
##   - cancer_sensitivity_comparison_v2.csv
##   - fig_14cancer_top50_v2.pdf
################################################################################

source("iwaba_functions.R")
source("iwaba_inference_helpers.R")

parse_arg_value <- function(args, prefix, default) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) == 0) {
    default
  } else {
    sub(paste0("^", prefix), "", hit[1])
  }
}

args <- commandArgs(trailingOnly = TRUE)
BOOT_B <- as.integer(parse_arg_value(args, "--B=", 999L))
ALPHA <- as.numeric(parse_arg_value(args, "--alpha=", 0.05))
BOOT_SEED <- as.integer(parse_arg_value(args, "--seed=", 20260319L))
MAX_PAIR_BOOT <- as.integer(parse_arg_value(args, "--max-pairs=", 10L))

DESC_VERSION <- "v2"
INF_VERSION <- "v1"

TOP50_LEVEL1_FILE <- sprintf("cancer_top50_level1_%s.csv", DESC_VERSION)
TOP50_LEVEL2_FILE <- sprintf("cancer_top50_level2_%s.csv", DESC_VERSION)
TOP50_LEVEL2_AGG_FILE <- sprintf("cancer_top50_level2_aggregate_%s.csv", DESC_VERSION)
TOP50_LEVEL3_FILE <- sprintf("cancer_top50_level3_%s.csv", DESC_VERSION)
TOP50_LEVEL4_FILE <- sprintf("cancer_top50_level4_%s.csv", DESC_VERSION)
TOP50_LEVEL5_FILE <- sprintf("cancer_top50_level5_%s.csv", DESC_VERSION)
TOP50_SUMMARY_FILE <- sprintf("cancer_top50_summary_%s.csv", DESC_VERSION)
TOP50_SCREENING_FILE <- sprintf("cancer_top50_screening_%s.csv", INF_VERSION)
TOP50_INFERENCE_FILE <- sprintf("cancer_top50_inference_%s.csv", INF_VERSION)
ALLGENES_SUMMARY_FILE <- sprintf("cancer_allgenes_summary_%s.csv", DESC_VERSION)
ALLGENES_SCREENING_FILE <- sprintf("cancer_allgenes_screening_%s.csv", INF_VERSION)
SENSITIVITY_FILE <- sprintf("cancer_sensitivity_comparison_%s.csv", DESC_VERSION)
FIG_FILE <- sprintf("fig_14cancer_top50_%s.pdf", DESC_VERSION)

CANCER_NAMES <- c(
  "Breast", "Prostate", "Lung", "Colorectal", "Lymphoma",
  "Bladder", "Melanoma", "Uterus", "Leukemia", "Renal",
  "Pancreas", "Ovary", "Mesothelioma", "CNS"
)

cat("================================================================\n")
cat("Loading 14-Cancer data (Ramaswamy et al., 2001, PNAS)\n")
cat("================================================================\n")

xtrain <- as.matrix(read.table("data/14cancer.xtrain"))
xtest <- as.matrix(read.table("data/14cancer.xtest"))
ytrain <- scan("data/14cancer.ytrain", what = integer(), quiet = TRUE)
ytest <- scan("data/14cancer.ytest", what = integer(), quiet = TRUE)

cat(sprintf("  Training: %d probes x %d samples\n", nrow(xtrain), ncol(xtrain)))
cat(sprintf("  Test:     %d probes x %d samples\n", nrow(xtest), ncol(xtest)))

X_all <- t(cbind(xtrain, xtest))
y_all <- c(ytrain, ytest)
cancer_labels <- factor(CANCER_NAMES[y_all], levels = CANCER_NAMES)

probe_ids <- paste0("Probe_", seq_len(ncol(X_all)))
colnames(X_all) <- probe_ids

cat(sprintf("  Combined: n = %d, p = %d, K = %d\n", nrow(X_all), ncol(X_all), nlevels(cancer_labels)))
cat(sprintf("  Bootstrap B = %d, alpha = %.3f, seed = %d, max pairwise bootstraps = %d\n",
            BOOT_B, ALPHA, BOOT_SEED, MAX_PAIR_BOOT))

cat("\n  Cancer type distribution:\n")
for (i in seq_along(CANCER_NAMES)) {
  cat(sprintf("    %2d. %-15s: n = %d\n", i, CANCER_NAMES[i], sum(cancer_labels == CANCER_NAMES[i])))
}

cat("\n================================================================\n")
cat("Preprocessing: variance filter + z-score normalization + ANOVA top-50\n")
cat("================================================================\n")

gene_var <- apply(X_all, 2, var)
X_all <- X_all[, gene_var > 1e-10, drop = FALSE]
cat(sprintf("  After QC: %d probes\n", ncol(X_all)))

X_scaled <- scale(X_all)
X_scaled <- as.matrix(X_scaled)
colnames(X_scaled) <- colnames(X_all)

grand_mean <- colMeans(X_scaled)
BSS <- numeric(ncol(X_scaled))
WSS <- numeric(ncol(X_scaled))

for (ct in levels(cancer_labels)) {
  mask <- cancer_labels == ct
  n_k <- sum(mask)
  group_mean <- colMeans(X_scaled[mask, , drop = FALSE])
  BSS <- BSS + n_k * (group_mean - grand_mean)^2
  WSS <- WSS + colSums(sweep(X_scaled[mask, , drop = FALSE], 2, group_mean)^2)
}

F_scores <- (BSS / (nlevels(cancer_labels) - 1)) / (WSS / (nrow(X_scaled) - nlevels(cancer_labels)))
F_scores[!is.finite(F_scores)] <- 0

top50_idx <- order(F_scores, decreasing = TRUE)[seq_len(50)]
top50_probes <- colnames(X_scaled)[top50_idx]
X_top50 <- X_scaled[, top50_idx, drop = FALSE]

cat(sprintf("  Top-50 F-score range: %.2f to %.2f\n",
            min(F_scores[top50_idx]), max(F_scores[top50_idx])))

build_level1_df <- function(result, X_matrix, varnames) {
  l1 <- result$level1
  eta <- l1$eta_results
  data.frame(
    variable = varnames,
    eta_B_waba = eta$eta_B_waba,
    eta_W_waba = eta$eta_W_waba,
    E_waba = eta$E_waba,
    eta_B_int = eta$eta_B_iwaba,
    eta_W_int = eta$eta_W_iwaba,
    E_int = eta$E_iwaba,
    delta_eta_B = eta$delta_eta,
    E_ratio = eta$E_ratio,
    class_waba = eta$class_waba,
    class_iwaba = eta$class_iwaba,
    V_C = l1$V_C,
    V_R = l1$V_R,
    V_W = l1$V_W,
    delta_E = iwaba_entity_dominance_contrast(result)[varnames],
    F_score = F_scores[match(varnames, colnames(X_scaled))],
    stringsAsFactors = FALSE
  )
}

build_level2_df <- function(result, varnames) {
  l2 <- result$level2
  pair_labels <- iwaba_pair_labels(l2$pairs, varnames)
  delta_A <- iwaba_pairwise_dominance_contrast(result, varnames = varnames)
  data.frame(
    pair_label = pair_labels,
    var1 = varnames[l2$pairs[1, ]],
    var2 = varnames[l2$pairs[2, ]],
    r_between_waba = l2$r_between_classical[cbind(l2$pairs[1, ], l2$pairs[2, ])],
    r_within = l2$r_within[cbind(l2$pairs[1, ], l2$pairs[2, ])],
    r_between_int = l2$r_between_int[cbind(l2$pairs[1, ], l2$pairs[2, ])],
    weighted_between_waba = l2$bw_waba,
    weighted_within_waba = l2$ww_waba,
    weighted_between_int = l2$bw_int,
    weighted_within_int = l2$ww_int,
    delta_A = unname(delta_A[pair_labels]),
    dominance_waba = ifelse(abs(l2$bw_waba) > abs(l2$ww_waba), "Between", "Within"),
    dominance_iwaba = ifelse(abs(l2$bw_int) > abs(l2$ww_int), "Between", "Within"),
    flipped = ifelse(abs(l2$bw_waba) > abs(l2$ww_waba), "Between", "Within") !=
      ifelse(abs(l2$bw_int) > abs(l2$ww_int), "Between", "Within"),
    stringsAsFactors = FALSE
  )
}

build_level2_aggregate_df <- function(result, level2_df) {
  data.frame(
    n_pairs = nrow(level2_df),
    mean_abs_r_between_waba = mean(abs(level2_df$r_between_waba)),
    mean_abs_r_between_int = mean(abs(level2_df$r_between_int)),
    mean_abs_weighted_between_waba = mean(abs(level2_df$weighted_between_waba)),
    mean_abs_weighted_within_waba = mean(abs(level2_df$weighted_within_waba)),
    mean_abs_weighted_between_int = mean(abs(level2_df$weighted_between_int)),
    mean_abs_weighted_within_int = mean(abs(level2_df$weighted_within_int)),
    proportion_between_dominant_waba = mean(level2_df$dominance_waba == "Between"),
    proportion_between_dominant_iwaba = mean(level2_df$dominance_iwaba == "Between"),
    consistency_rate_waba = result$level3$consistency_rate_waba,
    consistency_rate_iwaba = result$level3$consistency_rate_int,
    stringsAsFactors = FALSE
  )
}

build_level4_df <- function(result, varnames) {
  l4 <- result$level4
  data.frame(
    variable = varnames,
    pi_R = l4$pi_R,
    E_C_int = l4$E_C_int,
    E_R_int = l4$E_R_int,
    E_R_het_int = l4$E_R_het_int,
    V_R_het = l4$V_R_het,
    stringsAsFactors = FALSE
  )
}

build_level5_df <- function(result, varnames) {
  l2 <- result$level2
  l5 <- result$level5
  pair_labels <- iwaba_pair_labels(l2$pairs, varnames)
  data.frame(
    pair_label = pair_labels,
    var1 = varnames[l2$pairs[1, ]],
    var2 = varnames[l2$pairs[2, ]],
    r_R_billard = l5$r_R_billard,
    r_R_billard_defined = l5$r_R_billard_defined,
    r_R_het = l5$r_R_het,
    r_R_het_defined = l5$r_R_het_defined,
    stringsAsFactors = FALSE
  )
}

cat("\n================================================================\n")
cat("Top-50 branch: descriptive Levels I-V\n")
cat("================================================================\n")

result_top50 <- iwaba_full(X_top50, cancer_labels)
print_iwaba_summary(result_top50)

l1_top50 <- result_top50$level1
l4_top50 <- result_top50$level4
l5_top50 <- result_top50$level5

top50_level1_df <- build_level1_df(result_top50, X_top50, colnames(X_top50))
top50_level2_df <- build_level2_df(result_top50, colnames(X_top50))
top50_level2_agg_df <- build_level2_aggregate_df(result_top50, top50_level2_df)
top50_level3_df <- data.frame(
  consistency_rate_waba = result_top50$level3$consistency_rate_waba,
  consistency_rate_iwaba = result_top50$level3$consistency_rate_int,
  n_pairs = nrow(top50_level2_df),
  n_dominance_flips = sum(top50_level2_df$flipped),
  stringsAsFactors = FALSE
)
top50_level4_df <- build_level4_df(result_top50, colnames(X_top50))
top50_level5_df <- build_level5_df(result_top50, colnames(X_top50))

top50_summary_df <- data.frame(
  branch = "top50",
  n = nrow(X_top50),
  p = ncol(X_top50),
  K = nlevels(cancer_labels),
  G = l1_top50$info_gain,
  radius_prop = l1_top50$radius_prop,
  mean_eta_B_waba = mean(l1_top50$eta_results$eta_B_waba),
  mean_eta_B_int = mean(l1_top50$eta_results$eta_B_iwaba),
  mean_delta_eta_B = mean(l1_top50$eta_results$delta_eta),
  reclassified = sum(l1_top50$eta_results$class_waba != l1_top50$eta_results$class_iwaba),
  mean_pi_R = mean(l4_top50$pi_R),
  mean_E_R_het_int = mean(l4_top50$E_R_het_int),
  T_het = iwaba_global_t_het(l4_top50),
  mean_r_R_billard = l5_top50$mean_r_R_billard,
  mean_r_R_het_signed = l5_top50$mean_r_R_het_signed,
  mean_abs_r_between_waba = top50_level2_agg_df$mean_abs_r_between_waba,
  mean_abs_r_between_int = top50_level2_agg_df$mean_abs_r_between_int,
  mean_abs_weighted_between_waba = top50_level2_agg_df$mean_abs_weighted_between_waba,
  mean_abs_weighted_between_int = top50_level2_agg_df$mean_abs_weighted_between_int,
  mean_abs_weighted_within_waba = top50_level2_agg_df$mean_abs_weighted_within_waba,
  mean_abs_weighted_within_int = top50_level2_agg_df$mean_abs_weighted_within_int,
  consistency_rate_waba = result_top50$level3$consistency_rate_waba,
  consistency_rate_iwaba = result_top50$level3$consistency_rate_int,
  stringsAsFactors = FALSE
)

cat(sprintf("  Top-50 summary: G = %.3f, radius share = %.3f, mean signed r_R,het = %.3f\n",
            top50_summary_df$G, top50_summary_df$radius_prop, top50_summary_df$mean_r_R_het_signed))

cat("\nTop-50 inference: Brown-Forsythe screening + restricted bootstrap\n")
top50_screening_df <- iwaba_brown_forsythe_matrix(
  X_top50,
  cancer_labels,
  varnames = colnames(X_top50),
  alpha = ALPHA
)
top50_screening_df$F_score <- F_scores[match(top50_screening_df$variable, colnames(X_scaled))]
top50_screening_df$eta_B_int = l1_top50$eta_results$eta_B_iwaba[match(top50_screening_df$variable, colnames(X_top50))]
top50_screening_df$delta_eta_B = l1_top50$eta_results$delta_eta[match(top50_screening_df$variable, colnames(X_top50))]
top50_screening_df$pi_R = l4_top50$pi_R[match(top50_screening_df$variable, colnames(X_top50))]
top50_screening_df$E_R_het_int = l4_top50$E_R_het_int[match(top50_screening_df$variable, colnames(X_top50))]

screened_genes <- top50_screening_df$variable[top50_screening_df$reject]
candidate_pair_idx <- if (length(screened_genes) >= 2L) {
  which(top50_level2_df$var1 %in% screened_genes & top50_level2_df$var2 %in% screened_genes)
} else {
  integer(0)
}
if (length(candidate_pair_idx) == 0L) {
  candidate_pair_idx <- seq_len(nrow(top50_level2_df))
}
candidate_pair_idx <- candidate_pair_idx[order(abs(top50_level5_df$r_R_het[candidate_pair_idx]), decreasing = TRUE)]
selected_pair_idx <- head(candidate_pair_idx, min(MAX_PAIR_BOOT, length(candidate_pair_idx)))
selected_pair_labels <- top50_level2_df$pair_label[selected_pair_idx]

top50_boot_stat_fn <- function(Xb, gb) {
  result_b <- iwaba_full(Xb, gb)
  pair_labels_b <- iwaba_pair_labels(result_b$level2$pairs, colnames(X_top50))
  delta_A_b <- iwaba_pairwise_dominance_contrast(result_b, varnames = colnames(X_top50))
  names(delta_A_b) <- pair_labels_b
  r_R_het_b <- result_b$level5$r_R_het
  names(r_R_het_b) <- pair_labels_b

  stats::setNames(
    c(
      mean(result_b$level1$eta_results$eta_B_iwaba),
      mean(result_b$level4$pi_R),
      mean(result_b$level4$E_R_het_int),
      iwaba_global_t_het(result_b),
      delta_A_b[selected_pair_labels],
      r_R_het_b[selected_pair_labels]
    ),
    c(
      "aggregate::mean_eta_B_int",
      "aggregate::mean_pi_R",
      "aggregate::mean_E_R_het_int",
      "aggregate::T_het",
      paste0("delta_A::", selected_pair_labels),
      paste0("r_R_het::", selected_pair_labels)
    )
  )
}

top50_boot <- iwaba_within_group_bootstrap_vector(
  X = X_top50,
  g = cancer_labels,
  statistic_fn = top50_boot_stat_fn,
  B = BOOT_B,
  conf_level = 0.95,
  seed = BOOT_SEED,
  min_finite = max(50L, floor(BOOT_B * 0.50)),
  null_value = 0,
  progress = FALSE
)

top50_boot_ci <- top50_boot$ci
top50_name_parts <- strsplit(top50_boot_ci$stat_name, "::", fixed = TRUE)
top50_boot_ci$analysis_family <- vapply(top50_name_parts, `[`, character(1), 1)
top50_boot_ci$target_label <- vapply(top50_name_parts, `[`, character(1), 2)
top50_boot_ci$p_value <- NA_real_
top50_boot_ci$p_adj <- NA_real_
top50_boot_ci$reject <- NA
top50_boot_ci$notes <- ifelse(
  top50_boot_ci$analysis_family == "aggregate",
  "Within-group bootstrap CI for aggregate top-50 summary",
  ifelse(
    top50_boot_ci$analysis_family == "delta_A",
    "Restricted pairwise bootstrap CI for pairwise dominance contrast",
    "Restricted pairwise bootstrap CI for signed heterogeneity-based radius correlation"
  )
)

screening_summary_rows <- data.frame(
  analysis_family = "screening_summary",
  target_label = c("n_genes", "n_bf_raw_p_lt_alpha", "n_bh_reject", "prop_bh_reject", "selected_pair_count"),
  estimate = c(
    nrow(top50_screening_df),
    sum(top50_screening_df$p_value < ALPHA, na.rm = TRUE),
    sum(top50_screening_df$reject, na.rm = TRUE),
    mean(top50_screening_df$reject, na.rm = TRUE),
    length(selected_pair_labels)
  ),
  p_value = NA_real_,
  p_adj = NA_real_,
  reject = NA,
  ci_lower = NA_real_,
  ci_upper = NA_real_,
  ci_defined = NA,
  boot_defined_rate = NA_real_,
  reject_null = NA,
  notes = c(
    "Top-50 branch size",
    "Raw Brown-Forsythe discoveries",
    "BH-adjusted Brown-Forsythe discoveries",
    "BH-adjusted discovery proportion",
    "Restricted pairwise bootstrap subset size"
  ),
  stringsAsFactors = FALSE
)

top50_screening_export <- top50_screening_df
top50_screening_export$analysis_family <- "brown_forsythe"
top50_screening_export$target_label <- top50_screening_export$variable
top50_screening_export$estimate <- top50_screening_export$statistic
top50_screening_export$ci_lower <- NA_real_
top50_screening_export$ci_upper <- NA_real_
top50_screening_export$ci_defined <- NA
top50_screening_export$boot_defined_rate <- NA_real_
top50_screening_export$reject_null <- !is.na(top50_screening_export$p_value) & top50_screening_export$p_value < ALPHA
top50_screening_export$notes <- ifelse(
  top50_screening_export$reject,
  "BH-significant Brown-Forsythe heterogeneity screen",
  "Brown-Forsythe heterogeneity screen"
)

top50_inference_df <- rbind(
  top50_screening_export[, c(
    "analysis_family", "target_label", "estimate", "p_value", "p_adj", "reject",
    "ci_lower", "ci_upper", "ci_defined", "boot_defined_rate", "reject_null", "notes"
  )],
  screening_summary_rows,
  top50_boot_ci[, c(
    "analysis_family", "target_label", "estimate", "p_value", "p_adj", "reject",
    "ci_lower", "ci_upper", "ci_defined", "boot_defined_rate", "reject_null", "notes"
  )]
)

cat(sprintf("  Top-50 BF/BH discoveries: %d / %d\n",
            sum(top50_screening_df$reject), nrow(top50_screening_df)))
cat(sprintf("  Restricted pairwise bootstrap subset: %d pairs\n", length(selected_pair_labels)))

cat("\n================================================================\n")
cat("All-gene branch: descriptive Level I / IV + high-dimensional screening\n")
cat("================================================================\n")

result_all_l1 <- iwaba_level1(X_scaled, cancer_labels)
l4_all <- iwaba_level4(result_all_l1$gq, result_all_l1$V_C, result_all_l1$V_R, result_all_l1$V_W)
eta_all <- result_all_l1$eta_results

allgenes_screening_base <- iwaba_brown_forsythe_matrix(
  X_scaled,
  cancer_labels,
  varnames = colnames(X_scaled),
  alpha = ALPHA
)

allgenes_screening_df <- data.frame(
  variable = colnames(X_scaled),
  bf_statistic = allgenes_screening_base$statistic,
  bf_p_value = allgenes_screening_base$p_value,
  bf_p_adj = allgenes_screening_base$p_adj,
  bf_reject = allgenes_screening_base$reject,
  eta_B_waba = eta_all$eta_B_waba,
  eta_B_int = eta_all$eta_B_iwaba,
  delta_eta_B = eta_all$delta_eta,
  class_waba = eta_all$class_waba,
  class_iwaba = eta_all$class_iwaba,
  reclassified = eta_all$class_waba != eta_all$class_iwaba,
  pi_R = l4_all$pi_R,
  E_R_het_int = l4_all$E_R_het_int,
  V_R_het = l4_all$V_R_het,
  stringsAsFactors = FALSE
)

mean_radius_all <- rowMeans(result_all_l1$gq$radii)
allgenes_summary_df <- data.frame(
  branch = "all_genes",
  n = result_all_l1$n,
  p = result_all_l1$p,
  K = result_all_l1$K,
  G = result_all_l1$info_gain,
  radius_prop = result_all_l1$radius_prop,
  mean_eta_B_waba = mean(eta_all$eta_B_waba),
  mean_eta_B_int = mean(eta_all$eta_B_iwaba),
  mean_delta_eta_B = mean(eta_all$delta_eta),
  positive_delta_eta_count = sum(eta_all$delta_eta > 0),
  positive_delta_eta_prop = mean(eta_all$delta_eta > 0),
  reclassified_count = sum(eta_all$class_waba != eta_all$class_iwaba),
  reclassified_prop = mean(eta_all$class_waba != eta_all$class_iwaba),
  mean_pi_R = mean(l4_all$pi_R),
  mean_E_R_het_int = mean(l4_all$E_R_het_int),
  T_het = iwaba_global_t_het(l4_all),
  bf_raw_p_lt_alpha_count = sum(allgenes_screening_df$bf_p_value < ALPHA, na.rm = TRUE),
  bf_bh_reject_count = sum(allgenes_screening_df$bf_reject, na.rm = TRUE),
  bf_bh_reject_prop = mean(allgenes_screening_df$bf_reject, na.rm = TRUE),
  mean_within_group_sd_min = min(mean_radius_all),
  mean_within_group_sd_max = max(mean_radius_all),
  mean_within_group_sd_ratio = max(mean_radius_all) / min(mean_radius_all),
  stringsAsFactors = FALSE
)

sensitivity_df <- rbind(
  data.frame(
    branch = "top50",
    p = ncol(X_top50),
    G = top50_summary_df$G,
    radius_prop = top50_summary_df$radius_prop,
    mean_eta_B_waba = top50_summary_df$mean_eta_B_waba,
    mean_eta_B_int = top50_summary_df$mean_eta_B_int,
    mean_delta_eta_B = top50_summary_df$mean_delta_eta_B,
    reclassified_count = top50_summary_df$reclassified,
    reclassified_prop = top50_summary_df$reclassified / ncol(X_top50),
    mean_E_R_het_int = top50_summary_df$mean_E_R_het_int,
    stringsAsFactors = FALSE
  ),
  data.frame(
    branch = "all_genes",
    p = ncol(X_scaled),
    G = allgenes_summary_df$G,
    radius_prop = allgenes_summary_df$radius_prop,
    mean_eta_B_waba = allgenes_summary_df$mean_eta_B_waba,
    mean_eta_B_int = allgenes_summary_df$mean_eta_B_int,
    mean_delta_eta_B = allgenes_summary_df$mean_delta_eta_B,
    reclassified_count = allgenes_summary_df$reclassified_count,
    reclassified_prop = allgenes_summary_df$reclassified_prop,
    mean_E_R_het_int = allgenes_summary_df$mean_E_R_het_int,
    stringsAsFactors = FALSE
  )
)

write.csv(top50_level1_df, TOP50_LEVEL1_FILE, row.names = FALSE)
write.csv(top50_level2_df, TOP50_LEVEL2_FILE, row.names = FALSE)
write.csv(top50_level2_agg_df, TOP50_LEVEL2_AGG_FILE, row.names = FALSE)
write.csv(top50_level3_df, TOP50_LEVEL3_FILE, row.names = FALSE)
write.csv(top50_level4_df, TOP50_LEVEL4_FILE, row.names = FALSE)
write.csv(top50_level5_df, TOP50_LEVEL5_FILE, row.names = FALSE)
write.csv(top50_summary_df, TOP50_SUMMARY_FILE, row.names = FALSE)
write.csv(top50_screening_df, TOP50_SCREENING_FILE, row.names = FALSE)
write.csv(top50_inference_df, TOP50_INFERENCE_FILE, row.names = FALSE)
write.csv(allgenes_summary_df, ALLGENES_SUMMARY_FILE, row.names = FALSE)
write.csv(allgenes_screening_df, ALLGENES_SCREENING_FILE, row.names = FALSE)
write.csv(sensitivity_df, SENSITIVITY_FILE, row.names = FALSE)

cat("\nSaved outputs:\n")
cat(sprintf("  %s\n", TOP50_LEVEL1_FILE))
cat(sprintf("  %s\n", TOP50_LEVEL2_FILE))
cat(sprintf("  %s\n", TOP50_LEVEL2_AGG_FILE))
cat(sprintf("  %s\n", TOP50_LEVEL3_FILE))
cat(sprintf("  %s\n", TOP50_LEVEL4_FILE))
cat(sprintf("  %s\n", TOP50_LEVEL5_FILE))
cat(sprintf("  %s\n", TOP50_SUMMARY_FILE))
cat(sprintf("  %s\n", TOP50_SCREENING_FILE))
cat(sprintf("  %s\n", TOP50_INFERENCE_FILE))
cat(sprintf("  %s\n", ALLGENES_SUMMARY_FILE))
cat(sprintf("  %s\n", ALLGENES_SCREENING_FILE))
cat(sprintf("  %s\n", SENSITIVITY_FILE))

cat("\nGenerating top-50 figure...\n")
pdf(FIG_FILE, width = 14, height = 12)
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))

plot(l1_top50$eta_results$eta_B_waba, l1_top50$eta_results$eta_B_iwaba,
     xlab = expression(eta[B]^"(WABA)"),
     ylab = expression(eta[B]^"(I-WABA)"),
     main = "(a) Level I: Top-50 Entity Analysis",
     pch = 19, col = "steelblue", cex = 1.6,
     xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1, col = "red", lty = 2, lwd = 2.5)
grid(col = "gray80", lty = 1)

hist(l1_top50$eta_results$delta_eta, breaks = 15, col = "#DD8452",
     border = "darkgray",
     xlim = c(min(0, min(l1_top50$eta_results$delta_eta)), max(l1_top50$eta_results$delta_eta) * 1.05),
     main = expression("(b) Distribution of " * Delta * eta[B] * " (Top-50)"),
     xlab = expression(Delta * eta[B]),
     ylab = "Frequency")
abline(v = mean(l1_top50$eta_results$delta_eta), col = "blue", lwd = 2.5)
abline(v = 0, col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("Mean", "Zero"), col = c("blue", "red"), lty = c(1, 2), bty = "n")

mean_radius_top50 <- rowMeans(result_top50$gq$radii)
sorted_idx <- order(mean_radius_top50, decreasing = TRUE)
sorted_names <- names(sort(mean_radius_top50, decreasing = TRUE))
boxplot(t(result_top50$gq$radii[sorted_idx, , drop = FALSE]),
        las = 2, col = "#55A868",
        main = "(c) Within-Group SD by Cancer Type",
        ylab = "Within-group SD",
        cex.axis = 0.75,
        names = sorted_names,
        border = "darkgray", pch = 16)
grid(col = "gray80", lty = 1)

pi_R_sorted <- sort(l4_top50$pi_R, decreasing = TRUE)
barplot(pi_R_sorted,
        col = "#4C72B0",
        main = expression("(d) Level IV: Radius Share " * pi[R] * " (Top-50)"),
        ylab = expression(pi[R]),
        xlab = "Top-50 probes (sorted)",
        border = "darkgray",
        ylim = c(0, max(pi_R_sorted) * 1.1))
abline(h = mean(l4_top50$pi_R), col = "red", lty = 2, lwd = 2)
legend("topright", legend = "Mean", col = "red", lty = 2, bty = "n")
grid(col = "gray80", lty = 1)

dev.off()
cat(sprintf("  %s\n", FIG_FILE))

cat("\nKey summaries:\n")
cat(sprintf("  Top-50: G = %.3f, mean pi_R = %.3f, mean E_R,het = %.3f, mean signed r_R,het = %.3f\n",
            top50_summary_df$G, top50_summary_df$mean_pi_R,
            top50_summary_df$mean_E_R_het_int, top50_summary_df$mean_r_R_het_signed))
cat(sprintf("  All genes: G = %.3f, mean pi_R = %.3f, mean E_R,het = %.3f, T_het = %.3f\n",
            allgenes_summary_df$G, allgenes_summary_df$mean_pi_R,
            allgenes_summary_df$mean_E_R_het_int, allgenes_summary_df$T_het))
cat(sprintf("  All-gene BF/BH discoveries: %d / %d\n",
            allgenes_summary_df$bf_bh_reject_count, ncol(X_scaled)))

cat("\nDone.\n")
