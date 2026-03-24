################################################################################
## I-WABA Inferential Robustness Simulation
##
## Sensitivity block requested by reviewers:
##   - moderate skewness via standardized log-normal errors
##   - 10:1 group-size imbalance under Gaussian errors
##
## Targets:
##   - coverage and mean interval length for
##       delta_1^(E), delta_12^(A), r_R,het(X1, X2)
################################################################################

source("iwaba_functions.R")
source("iwaba_inference_helpers.R")

OUTPUT_VERSION <- "v1"
RESULTS_DIR <- "results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

DETAIL_RDS <- file.path(RESULTS_DIR, sprintf("simulation_inference_robustness_results_%s.rds", OUTPUT_VERSION))
DETAIL_CSV <- file.path(RESULTS_DIR, sprintf("simulation_inference_robustness_results_%s.csv", OUTPUT_VERSION))
SUMMARY_RDS <- file.path(RESULTS_DIR, sprintf("simulation_inference_robustness_summary_%s.rds", OUTPUT_VERSION))
SUMMARY_CSV <- file.path(RESULTS_DIR, sprintf("simulation_inference_robustness_summary_%s.csv", OUTPUT_VERSION))

DEFAULT_R <- 100L
DEFAULT_B <- 199L
ALPHA <- 0.05
K_FIXED <- 5L
P_FIXED <- 5L
N_VALUES <- c(100L, 500L)
SETTINGS <- c("skewed_lognormal_balanced", "gaussian_unbalanced")

parse_arg_value <- function(args, prefix, default) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) == 0L) default else sub(paste0("^", prefix), "", hit[1])
}

args <- commandArgs(trailingOnly = TRUE)
MC_REPS <- as.integer(parse_arg_value(args, "--R=", DEFAULT_R))
BOOT_REPS <- as.integer(parse_arg_value(args, "--B=", DEFAULT_B))

if (is.na(MC_REPS) || MC_REPS <= 0L) stop("MC replications must be a positive integer.")
if (is.na(BOOT_REPS) || BOOT_REPS <= 0L) stop("Bootstrap resamples must be a positive integer.")

mean_pattern <- c(-1.00, -0.50, 0.00, 0.50, 1.00)
mean_scales <- c(1.40, 1.60, 1.45, 1.80, 1.35)

make_aligned_template <- function(weights) {
  means <- outer(mean_pattern, mean_scales)
  colnames(means) <- paste0("X", seq_len(P_FIXED))
  rownames(means) <- paste0("G", seq_len(K_FIXED))

  pattern_x1 <- c(-0.36, -0.18, 0.00, 0.18, 0.36)
  pattern_x2 <- c(-0.24, -0.18, 0.12, 0.12, 0.18)
  sds <- cbind(
    0.95 + pattern_x1,
    1.05 + pattern_x2,
    1.10 + 0.80 * pattern_x1,
    1.00 + 0.60 * pattern_x2,
    1.15 + 0.70 * pattern_x1
  )
  colnames(sds) <- paste0("X", seq_len(P_FIXED))
  rownames(sds) <- paste0("G", seq_len(K_FIXED))

  list(means = means, sds = sds, weights = weights)
}

weighted_center <- function(M, weights) {
  col_means <- colSums(sweep(M, 1, weights, "*"))
  sweep(M, 2, col_means)
}

weighted_corr_pair <- function(x, y, weights, eps = 1e-15) {
  x_centered <- x - sum(weights * x)
  y_centered <- y - sum(weights * y)
  denom <- sqrt(sum(weights * x_centered^2) * sum(weights * y_centered^2))
  if (denom <= eps) NA_real_ else max(min(sum(weights * x_centered * y_centered) / denom, 1), -1)
}

true_population_targets <- function(template) {
  weights <- template$weights
  means <- template$means
  sds <- template$sds
  true_rhet <- weighted_corr_pair(sds[, 1], sds[, 2], weights)

  mean_centered <- weighted_center(means, weights)
  V_C <- colSums(sweep(mean_centered^2, 1, weights, "*"))
  V_R <- (1 / 3) * colSums(sweep(sds^2, 1, weights, "*"))
  V_W <- colSums(sweep(sds^2, 1, weights, "*"))

  Cov_between_centers <- matrix(0, nrow = P_FIXED, ncol = P_FIXED)
  Cov_between_radii <- matrix(0, nrow = P_FIXED, ncol = P_FIXED)
  for (k in seq_len(K_FIXED)) {
    Cov_between_centers <- Cov_between_centers + weights[k] * tcrossprod(mean_centered[k, ])
    Cov_between_radii <- Cov_between_radii + weights[k] * tcrossprod(sds[k, ])
  }
  Cov_between_int <- Cov_between_centers + Cov_between_radii / 3
  Cov_within <- diag(V_W, nrow = P_FIXED, ncol = P_FIXED)
  Cov_dec <- Cov_between_int + Cov_within

  eta_B_int <- sqrt(safe_div(diag(Cov_between_int), diag(Cov_dec)))
  eta_W_int <- sqrt(safe_div(diag(Cov_within), diag(Cov_dec)))
  r_between_int <- cov_to_corr(Cov_between_int)
  r_within <- cov_to_corr(Cov_within)

  D_B_12 <- eta_B_int[1] * eta_B_int[2] * r_between_int[1, 2]
  D_W_12 <- eta_W_int[1] * eta_W_int[2] * r_within[1, 2]

  list(
    delta_E_1 = V_C[1] + V_R[1] - V_W[1],
    delta_A_12 = D_B_12^2 - D_W_12^2,
    r_R_het_12 = true_rhet
  )
}

balanced_sizes <- function(n) {
  rep.int(n %/% K_FIXED, K_FIXED)
}

unbalanced_sizes <- function(n) {
  ratios <- c(1, 2, 3, 5, 10)
  raw_sizes <- n * ratios / sum(ratios)
  sizes <- floor(raw_sizes)
  remainder <- n - sum(sizes)
  if (remainder > 0L) {
    add_idx <- order(raw_sizes - sizes, decreasing = TRUE)[seq_len(remainder)]
    sizes[add_idx] <- sizes[add_idx] + 1L
  }
  sizes
}

standardized_lognormal <- function(n, sdlog = 1) {
  pop_mean <- exp(sdlog^2 / 2)
  pop_sd <- sqrt((exp(sdlog^2) - 1) * exp(sdlog^2))
  (rlnorm(n, meanlog = 0, sdlog = sdlog) - pop_mean) / pop_sd
}

simulate_setting <- function(setting, n, seed) {
  set.seed(seed)

  group_sizes <- if (setting == "gaussian_unbalanced") unbalanced_sizes(n) else balanced_sizes(n)
  weights <- group_sizes / sum(group_sizes)
  template <- make_aligned_template(weights)

  X <- matrix(NA_real_, nrow = sum(group_sizes), ncol = P_FIXED)
  colnames(X) <- paste0("X", seq_len(P_FIXED))
  Y <- factor(rep(seq_len(K_FIXED), times = group_sizes), levels = seq_len(K_FIXED))

  row_start <- 1L
  for (k in seq_len(K_FIXED)) {
    nk <- group_sizes[k]
    row_end <- row_start + nk - 1L
    for (j in seq_len(P_FIXED)) {
      z <- if (setting == "skewed_lognormal_balanced") {
        standardized_lognormal(nk, sdlog = 1)
      } else {
        rnorm(nk)
      }
      X[row_start:row_end, j] <- template$means[k, j] + template$sds[k, j] * z
    }
    row_start <- row_end + 1L
  }

  list(
    X = X,
    Y = Y,
    group_sizes = group_sizes,
    truth = true_population_targets(template)
  )
}

extract_targets <- function(analysis) {
  pair_idx <- which(analysis$level2$pairs[1, ] == 1L & analysis$level2$pairs[2, ] == 2L)
  if (length(pair_idx) != 1L) stop("Could not locate pair (1,2).")

  stats::setNames(
    c(
      analysis$level1$V_C[1] + analysis$level1$V_R[1] - analysis$level1$V_W[1],
      analysis$level2$bw_int[pair_idx]^2 - analysis$level2$ww_int[pair_idx]^2,
      analysis$level5$r_R_het[pair_idx]
    ),
    c("delta_E_1", "delta_A_12", "r_R_het_12")
  )
}

make_seed <- function(setting_index, n, rep) {
  setting_index * 1000000L + n * 1000L + rep
}

rows <- vector("list", length(SETTINGS) * length(N_VALUES) * MC_REPS)
counter <- 1L

for (setting_index in seq_along(SETTINGS)) {
  setting <- SETTINGS[setting_index]
  for (n in N_VALUES) {
    for (rep in seq_len(MC_REPS)) {
      data_seed <- make_seed(setting_index, n, rep)
      dat <- simulate_setting(setting, n, data_seed)
      analysis <- iwaba_full(dat$X, dat$Y)

      stat_fn <- function(Xb, gb) {
        extract_targets(iwaba_full(Xb, gb))
      }

      boot <- iwaba_within_group_bootstrap_vector(
        X = dat$X,
        g = dat$Y,
        statistic_fn = stat_fn,
        B = BOOT_REPS,
        conf_level = 1 - ALPHA,
        seed = data_seed + 500000L,
        min_finite = 50L,
        null_value = 0,
        progress = FALSE
      )

      ci <- boot$ci
      truth <- dat$truth
      truth_vec <- c(truth$delta_E_1, truth$delta_A_12, truth$r_R_het_12)

      row <- data.frame(
        setting = setting,
        n = n,
        rep = rep,
        B = BOOT_REPS,
        delta_E_truth = truth_vec[1],
        delta_A_truth = truth_vec[2],
        r_R_het_truth = truth_vec[3],
        stringsAsFactors = FALSE
      )

      for (j in seq_len(nrow(ci))) {
        stat_name <- ci$stat_name[j]
        prefix <- switch(
          stat_name,
          "delta_E_1" = "delta_E",
          "delta_A_12" = "delta_A",
          "r_R_het_12" = "r_R_het",
          stat_name
        )

        row[[paste0(prefix, "_estimate")]] <- ci$estimate[j]
        row[[paste0(prefix, "_lower")]] <- ci$ci_lower[j]
        row[[paste0(prefix, "_upper")]] <- ci$ci_upper[j]
        row[[paste0(prefix, "_defined")]] <- ci$ci_defined[j]
        row[[paste0(prefix, "_length")]] <- ci$ci_upper[j] - ci$ci_lower[j]
        row[[paste0(prefix, "_coverage")]] <- isTRUE(ci$ci_defined[j]) &&
          is.finite(truth_vec[j]) &&
          truth_vec[j] >= ci$ci_lower[j] &&
          truth_vec[j] <= ci$ci_upper[j]
      }

      row$group_sizes <- paste(dat$group_sizes, collapse = ",")
      rows[[counter]] <- row
      counter <- counter + 1L
    }
  }
}

results <- do.call(rbind, rows)

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

summary_rows <- vector("list", length(SETTINGS) * length(N_VALUES))
counter <- 1L
for (setting in SETTINGS) {
  for (n in N_VALUES) {
    subset_rows <- results[results$setting == setting & results$n == n, , drop = FALSE]
    summary_rows[[counter]] <- data.frame(
      setting = setting,
      n = n,
      mc_reps = nrow(subset_rows),
      B = BOOT_REPS,
      delta_E_coverage = mean_or_na(subset_rows$delta_E_coverage),
      delta_E_length = mean_or_na(subset_rows$delta_E_length),
      delta_A_coverage = mean_or_na(subset_rows$delta_A_coverage),
      delta_A_length = mean_or_na(subset_rows$delta_A_length),
      r_R_het_coverage = mean_or_na(subset_rows$r_R_het_coverage),
      r_R_het_length = mean_or_na(subset_rows$r_R_het_length),
      r_R_het_defined_rate = mean_or_na(as.numeric(subset_rows$r_R_het_defined)),
      example_group_sizes = subset_rows$group_sizes[1],
      stringsAsFactors = FALSE
    )
    counter <- counter + 1L
  }
}

summary_table <- do.call(rbind, summary_rows)

write.csv(results, DETAIL_CSV, row.names = FALSE)
saveRDS(results, DETAIL_RDS)
write.csv(summary_table, SUMMARY_CSV, row.names = FALSE)
saveRDS(summary_table, SUMMARY_RDS)

cat("Robustness simulation complete.\n")
print(summary_table, row.names = FALSE)
