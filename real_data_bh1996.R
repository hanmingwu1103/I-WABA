################################################################################
## I-WABA Real Data Analysis: Bliese-Halverson (bh1996) Dataset
## Revised descriptive + inferential workflow
##
## Outputs:
##   - bh1996_level1_v2.csv
##   - bh1996_level2_v2.csv
##   - bh1996_level3_v2.csv
##   - bh1996_level4_v2.csv
##   - bh1996_level5_v2.csv
##   - bh1996_summary_v2.csv
##   - bh1996_inference_v1.csv
##   - fig_bh1996_analysis_v2.pdf
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

DESC_VERSION <- "v2"
INF_VERSION <- "v1"

LEVEL1_FILE <- sprintf("bh1996_level1_%s.csv", DESC_VERSION)
LEVEL2_FILE <- sprintf("bh1996_level2_%s.csv", DESC_VERSION)
LEVEL3_FILE <- sprintf("bh1996_level3_%s.csv", DESC_VERSION)
LEVEL4_FILE <- sprintf("bh1996_level4_%s.csv", DESC_VERSION)
LEVEL5_FILE <- sprintf("bh1996_level5_%s.csv", DESC_VERSION)
SUMMARY_FILE <- sprintf("bh1996_summary_%s.csv", DESC_VERSION)
INFERENCE_FILE <- sprintf("bh1996_inference_%s.csv", INF_VERSION)
FIG_FILE <- sprintf("fig_bh1996_analysis_%s.pdf", DESC_VERSION)

cat("Loading bh1996 data...\n")

if (require(multilevel, quietly = TRUE)) {
  data(bh1996)
  df <- bh1996
  cat("  Loaded from multilevel R package\n")
} else {
  df <- read.csv("data/bh1996_rdata.csv", row.names = 1)
  cat("  Loaded from data/bh1996_rdata.csv\n")
}

variables <- c("COHES", "LEAD", "WBEING", "HRS")
group_col <- "GRP"

X <- as.matrix(df[, variables])
g <- as.factor(df[[group_col]])

cat(sprintf("  n = %d, K = %d, p = %d\n", nrow(X), nlevels(g), ncol(X)))
cat(sprintf("  Bootstrap B = %d, alpha = %.3f, seed = %d\n", BOOT_B, ALPHA, BOOT_SEED))

result <- iwaba_full(X, g)
print_iwaba_summary(result)

l1 <- result$level1
l2 <- result$level2
l3 <- result$level3
l4 <- result$level4
l5 <- result$level5
eta <- l1$eta_results
pair_labels <- iwaba_pair_labels(l2$pairs, variables)
delta_E <- iwaba_entity_dominance_contrast(result)
delta_A <- iwaba_pairwise_dominance_contrast(result, varnames = variables)

icc_values <- vapply(variables, function(v) {
  waba <- classical_waba(X[, v, drop = FALSE], g)
  unname(diag(waba$BSS) / diag(waba$TSS))
}, numeric(1))

level1_df <- data.frame(
  variable = variables,
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
  delta_E = unname(delta_E[variables]),
  ICC1 = unname(icc_values[variables]),
  stringsAsFactors = FALSE
)

level2_df <- data.frame(
  pair_label = pair_labels,
  var1 = variables[l2$pairs[1, ]],
  var2 = variables[l2$pairs[2, ]],
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

level3_df <- data.frame(
  consistency_rate_waba = l3$consistency_rate_waba,
  consistency_rate_iwaba = l3$consistency_rate_int,
  n_pairs = ncol(l2$pairs),
  n_dominance_flips = sum(level2_df$flipped),
  stringsAsFactors = FALSE
)

level4_df <- data.frame(
  variable = variables,
  pi_R = l4$pi_R,
  E_C_int = l4$E_C_int,
  E_R_int = l4$E_R_int,
  E_R_het_int = l4$E_R_het_int,
  V_R_het = l4$V_R_het,
  stringsAsFactors = FALSE
)

level5_df <- data.frame(
  pair_label = pair_labels,
  var1 = variables[l2$pairs[1, ]],
  var2 = variables[l2$pairs[2, ]],
  r_R_billard = l5$r_R_billard,
  r_R_billard_defined = l5$r_R_billard_defined,
  r_R_het = l5$r_R_het,
  r_R_het_defined = l5$r_R_het_defined,
  stringsAsFactors = FALSE
)

summary_df <- data.frame(
  dataset = "bh1996",
  n = l1$n,
  K = l1$K,
  p = l1$p,
  G = l1$info_gain,
  radius_prop = l1$radius_prop,
  mean_eta_B_waba = mean(eta$eta_B_waba),
  mean_eta_B_iwaba = mean(eta$eta_B_iwaba),
  mean_delta_eta_B = mean(eta$delta_eta),
  reclassified = sum(eta$class_waba != eta$class_iwaba),
  mean_pi_R = mean(l4$pi_R),
  mean_E_R_het_int = mean(l4$E_R_het_int),
  mean_r_R_billard = l5$mean_r_R_billard,
  mean_r_R_het_signed = l5$mean_r_R_het_signed,
  T_het = iwaba_global_t_het(l4),
  level3_consistency_waba = l3$consistency_rate_waba,
  level3_consistency_iwaba = l3$consistency_rate_int,
  level2_dominance_flips = sum(level2_df$flipped),
  stringsAsFactors = FALSE
)

bf_df <- iwaba_brown_forsythe_matrix(X, g, varnames = variables, alpha = ALPHA)
bf_df$analysis_family <- "brown_forsythe"
bf_df$target_label <- bf_df$variable
bf_df$estimate <- bf_df$statistic
bf_df$ci_lower <- NA_real_
bf_df$ci_upper <- NA_real_
bf_df$ci_defined <- NA
bf_df$boot_defined_rate <- NA_real_
bf_df$reject_null <- !is.na(bf_df$p_value) & bf_df$p_value < ALPHA
bf_df$notes <- "H0: equal within-group variances"

bh_boot_stat_fn <- function(Xb, gb) {
  result_b <- iwaba_full(Xb, gb)
  l4_b <- result_b$level4
  delta_E_b <- iwaba_entity_dominance_contrast(result_b)
  delta_A_b <- iwaba_pairwise_dominance_contrast(result_b, varnames = variables)
  names(delta_A_b) <- pair_labels
  r_R_het_b <- result_b$level5$r_R_het
  names(r_R_het_b) <- pair_labels

  stats::setNames(
    c(
      delta_E_b[variables],
      l4_b$pi_R,
      l4_b$E_R_het_int,
      delta_A_b[pair_labels],
      r_R_het_b[pair_labels]
    ),
    c(
      paste0("delta_E::", variables),
      paste0("pi_R::", variables),
      paste0("E_R_het::", variables),
      paste0("delta_A::", pair_labels),
      paste0("r_R_het::", pair_labels)
    )
  )
}

cat("\nRunning bh1996 bootstrap inference...\n")
bh_boot <- iwaba_within_group_bootstrap_vector(
  X = X,
  g = g,
  statistic_fn = bh_boot_stat_fn,
  B = BOOT_B,
  conf_level = 0.95,
  seed = BOOT_SEED,
  min_finite = max(50L, floor(BOOT_B * 0.50)),
  null_value = 0,
  progress = FALSE
)

boot_ci_df <- bh_boot$ci
name_parts <- strsplit(boot_ci_df$stat_name, "::", fixed = TRUE)
boot_ci_df$analysis_family <- vapply(name_parts, `[`, character(1), 1)
boot_ci_df$target_label <- vapply(name_parts, `[`, character(1), 2)
boot_ci_df$p_value <- NA_real_
boot_ci_df$p_adj <- NA_real_
boot_ci_df$reject <- NA
boot_ci_df$notes <- ifelse(
  boot_ci_df$analysis_family == "delta_E",
  "Bootstrap CI for entity-level dominance contrast",
  ifelse(
    boot_ci_df$analysis_family == "pi_R",
    "Bootstrap CI for radius share",
    ifelse(
      boot_ci_df$analysis_family == "E_R_het",
      "Bootstrap CI for heterogeneity-focused dispersion index",
      ifelse(
        boot_ci_df$analysis_family == "delta_A",
        "Bootstrap CI for pairwise dominance contrast",
        ifelse(
          boot_ci_df$analysis_family == "r_R_het",
          "Bootstrap CI for signed heterogeneity-based radius correlation",
          ""
        )
      )
    )
  )
)

inference_df <- rbind(
  bf_df[, c(
    "analysis_family", "target_label", "estimate", "p_value", "p_adj", "reject",
    "ci_lower", "ci_upper", "ci_defined", "boot_defined_rate", "reject_null", "notes"
  )],
  boot_ci_df[, c(
    "analysis_family", "target_label", "estimate", "p_value", "p_adj", "reject",
    "ci_lower", "ci_upper", "ci_defined", "boot_defined_rate", "reject_null", "notes"
  )]
)

write.csv(level1_df, LEVEL1_FILE, row.names = FALSE)
write.csv(level2_df, LEVEL2_FILE, row.names = FALSE)
write.csv(level3_df, LEVEL3_FILE, row.names = FALSE)
write.csv(level4_df, LEVEL4_FILE, row.names = FALSE)
write.csv(level5_df, LEVEL5_FILE, row.names = FALSE)
write.csv(summary_df, SUMMARY_FILE, row.names = FALSE)
write.csv(inference_df, INFERENCE_FILE, row.names = FALSE)

cat("\nSaved descriptive outputs:\n")
cat(sprintf("  %s\n", LEVEL1_FILE))
cat(sprintf("  %s\n", LEVEL2_FILE))
cat(sprintf("  %s\n", LEVEL3_FILE))
cat(sprintf("  %s\n", LEVEL4_FILE))
cat(sprintf("  %s\n", LEVEL5_FILE))
cat(sprintf("  %s\n", SUMMARY_FILE))
cat(sprintf("  %s\n", INFERENCE_FILE))

cat("\nGenerating figure...\n")
pdf(FIG_FILE, width = 14, height = 10)
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))

plot(eta$eta_B_waba, eta$eta_B_iwaba,
     xlab = expression(eta[B]^"(WABA)"),
     ylab = expression(eta[B]^"(I-WABA)"),
     main = "(a) Level I Entity Analysis",
     pch = 19, col = "steelblue", cex = 2.5,
     xlim = c(0.15, 0.5), ylim = c(0.15, 0.65))
abline(0, 1, col = "red", lty = 2, lwd = 2.5)
text(eta$eta_B_waba, eta$eta_B_iwaba, labels = eta$variable,
     pos = 3, cex = 1.1)
grid(col = "gray80", lty = 1)

boxplot(result$gq$radii,
        main = "(b) Within-Group SD Distributions",
        names = variables, col = "#55A868",
        ylab = "Within-group SD",
        border = "darkgray", pch = 16, cex = 0.8)
grid(col = "gray80", lty = 1)

barplot(rbind(l4$E_C_int, l4$E_R_het_int), beside = TRUE,
        names.arg = variables,
        col = c("steelblue", "coral"),
        main = "(c) Level IV: Center vs Heterogeneity",
        ylab = "E-ratio",
        ylim = c(0, max(l4$E_C_int, l4$E_R_het_int) * 1.15),
        border = "darkgray")
legend("topleft",
       legend = c(expression(E[C]^"(Int)"), expression(E[R*","*het]^"(Int)")),
       fill = c("steelblue", "coral"), bty = "n")
grid(col = "gray80", lty = 1)

level5_mat <- rbind(l5$r_R_billard, l5$r_R_het)
y_min <- min(0, level5_mat, na.rm = TRUE) * 1.15
y_max <- max(level5_mat, na.rm = TRUE) * 1.15
barplot(level5_mat, beside = TRUE,
        names.arg = gsub("--", "\n", pair_labels, fixed = TRUE),
        las = 1,
        col = c("steelblue", "coral"),
        main = "(d) Level V: Radius Associations",
        ylab = "Association",
        ylim = c(y_min, y_max),
        border = "darkgray")
abline(h = 0, col = "gray40", lwd = 1)
legend("topleft",
       legend = c(expression(r[R]^"(Billard)"), expression(r[R*","*het])),
       fill = c("steelblue", "coral"), bty = "n")
grid(col = "gray80", lty = 1)

dev.off()
cat(sprintf("  %s\n", FIG_FILE))

cat("\nKey descriptive summaries:\n")
cat(sprintf("  G = %.3f, radius share = %.3f, mean E_R,het = %.3f\n",
            summary_df$G, summary_df$radius_prop, summary_df$mean_E_R_het_int))
cat(sprintf("  Level II flips = %d / %d\n",
            summary_df$level2_dominance_flips, level3_df$n_pairs))
cat(sprintf("  Mean r_R^(Billard) = %.3f, mean signed r_R,het = %.3f\n",
            summary_df$mean_r_R_billard, summary_df$mean_r_R_het_signed))

cat("\nDone.\n")
