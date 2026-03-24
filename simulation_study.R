################################################################################
## I-WABA Simulation Study (v4)
## Title: Interval-based Within-And-Between Analysis
## Author: Han-Ming Wu (National Chengchi University)
##
## Revised Section 4 design:
##   - five scenarios, including complete_null
##   - matched-template p-sensitivity design
##   - explicit equivocal band thresholds
##   - K = 2 sign-only / structurally degenerate handling for r_R_het
##   - targeted comparator experiment for spread_only and complete_null
##   - versioned results/ outputs for manuscript use
################################################################################

# ============================================================================
# 0. Setup
# ============================================================================
if (!require("MASS"))      install.packages("MASS")
if (!require("ggplot2"))   install.packages("ggplot2")
if (!require("dplyr"))     install.packages("dplyr")
if (!require("tidyr"))     install.packages("tidyr")
if (!require("gridExtra")) install.packages("gridExtra")

library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

source("iwaba_functions.R")
source("iwaba_inference_helpers.R")

RESULTS_DIR <- "results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

OUTPUT_VERSION <- "v5"
METADATA_VERSION <- "v2"
COMPLETE_NULL_VERSION <- "v2"
COMPARATOR_VERSION <- "v2"
P_SENSITIVITY_VERSION <- "v3"
HEATMAP_VERSION <- "v3"

RESULTS_RDS <- file.path(RESULTS_DIR, sprintf("simulation_results_%s.rds", OUTPUT_VERSION))
RESULTS_CSV <- file.path(RESULTS_DIR, sprintf("simulation_results_%s.csv", OUTPUT_VERSION))
METADATA_RDS <- file.path(RESULTS_DIR, sprintf("simulation_metadata_%s.rds", METADATA_VERSION))
METADATA_CSV <- file.path(RESULTS_DIR, sprintf("simulation_metadata_%s.csv", METADATA_VERSION))
FOCAL_SUMMARY_RDS <- file.path(RESULTS_DIR, sprintf("simulation_focal_summary_%s.rds", OUTPUT_VERSION))
FOCAL_SUMMARY_CSV <- file.path(RESULTS_DIR, sprintf("simulation_focal_summary_%s.csv", OUTPUT_VERSION))
P_SENSITIVITY_RDS <- file.path(RESULTS_DIR, sprintf("simulation_p_sensitivity_%s.rds", P_SENSITIVITY_VERSION))
P_SENSITIVITY_CSV <- file.path(RESULTS_DIR, sprintf("simulation_p_sensitivity_%s.csv", P_SENSITIVITY_VERSION))
HEATMAP_RDS <- file.path(RESULTS_DIR, sprintf("simulation_heatmap_%s.rds", HEATMAP_VERSION))
HEATMAP_CSV <- file.path(RESULTS_DIR, sprintf("simulation_heatmap_%s.csv", HEATMAP_VERSION))
COMPLETE_NULL_RDS <- file.path(RESULTS_DIR, sprintf("simulation_complete_null_%s.rds", COMPLETE_NULL_VERSION))
COMPLETE_NULL_CSV <- file.path(RESULTS_DIR, sprintf("simulation_complete_null_%s.csv", COMPLETE_NULL_VERSION))
COMPARATOR_RDS <- file.path(RESULTS_DIR, sprintf("simulation_comparator_%s.rds", COMPARATOR_VERSION))
COMPARATOR_CSV <- file.path(RESULTS_DIR, sprintf("simulation_comparator_%s.csv", COMPARATOR_VERSION))

FIG_INFO_GAIN_FILE <- file.path(RESULTS_DIR, sprintf("fig_info_gain_vs_K_%s.pdf", OUTPUT_VERSION))
FIG_ETA_IMPROVEMENT_FILE <- file.path(RESULTS_DIR, sprintf("fig_eta_improvement_%s.pdf", OUTPUT_VERSION))
FIG_DISPERSION_FILE <- file.path(RESULTS_DIR, sprintf("fig_dispersion_diagnostics_%s.pdf", OUTPUT_VERSION))
FIG_HEATMAP_FILE <- file.path(RESULTS_DIR, sprintf("fig_heatmap_info_gain_%s.pdf", OUTPUT_VERSION))
FIG_RADIUS_FILE <- file.path(RESULTS_DIR, sprintf("fig_radius_association_%s.pdf", OUTPUT_VERSION))

SIM_N_REPS <- 100
SIM_K_VALUES <- c(2, 3, 5, 10)
SIM_N_VALUES <- c(100, 200, 500, 1000)
SIM_P_VALUES <- c(2, 5, 10, 50)
SIM_SCENARIOS <- c("heteroscedastic", "spread_only", "homoscedastic", "nonlinear_var", "complete_null")

FOCAL_N <- 500
FOCAL_P <- 10
P_SENSITIVITY_K <- 5
P_SENSITIVITY_N <- 500
COMPARATOR_SCENARIOS <- c("spread_only", "complete_null")
COMPARATOR_ALPHA <- 0.05

EQUIVOCAL_WITHIN_DEG <- 37.5
EQUIVOCAL_BETWEEN_DEG <- 52.5
EQUIVOCAL_WITHIN_THRESHOLD <- tan(EQUIVOCAL_WITHIN_DEG * pi / 180)
EQUIVOCAL_BETWEEN_THRESHOLD <- tan(EQUIVOCAL_BETWEEN_DEG * pi / 180)

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

sd_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) NA_real_ else stats::sd(x)
}

mcse_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    NA_real_
  } else if (length(x) == 1) {
    0
  } else {
    stats::sd(x) / sqrt(length(x))
  }
}

min_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else min(x)
}

format_num_or_dash <- function(x, digits = 3) {
  if (is.na(x)) {
    "--"
  } else {
    sprintf(paste0("%.", digits, "f"), x)
  }
}

classify_e_explicit <- function(E) {
  ifelse(E > EQUIVOCAL_BETWEEN_THRESHOLD, "Between",
         ifelse(E < EQUIVOCAL_WITHIN_THRESHOLD, "Within", "Equivocal"))
}

summarize_level5_mode <- function(x) {
  modes <- unique(x)
  if (length(modes) == 1) modes else "mixed"
}

format_level5_note <- function(mode, prop_pos, prop_neg, prop_undef) {
  if (mode == "k2_sign_only") {
    sprintf("K=2 sign-only: P(+)=%.2f, P(-)=%.2f, undef=%.2f",
            prop_pos, prop_neg, prop_undef)
  } else if (mode == "undefined") {
    sprintf("Undefined: centered radii vanish (undef=%.2f)", prop_undef)
  } else {
    ""
  }
}

subset_sim_data <- function(dat, p) {
  list(X = dat$X[, seq_len(p), drop = FALSE], Y = dat$Y)
}

make_sim_seed <- function(scenario_index, K, n, rep) {
  scenario_index * 1000000 + K * 10000 + n * 10 + rep
}

simulation_metadata_df <- function() {
  data.frame(
    output_version = OUTPUT_VERSION,
    metadata_version = METADATA_VERSION,
    complete_null_version = COMPLETE_NULL_VERSION,
    comparator_version = COMPARATOR_VERSION,
    p_sensitivity_version = P_SENSITIVITY_VERSION,
    heatmap_version = HEATMAP_VERSION,
    n_reps = SIM_N_REPS,
    K_values = paste(SIM_K_VALUES, collapse = ","),
    n_values = paste(SIM_N_VALUES, collapse = ","),
    p_values = paste(SIM_P_VALUES, collapse = ","),
    scenarios = paste(SIM_SCENARIOS, collapse = ","),
    focal_n = FOCAL_N,
    focal_p = FOCAL_P,
    p_sensitivity_K = P_SENSITIVITY_K,
    p_sensitivity_n = P_SENSITIVITY_N,
    comparator_scenarios = paste(COMPARATOR_SCENARIOS, collapse = ","),
    comparator_alpha = COMPARATOR_ALPHA,
    equivocal_within_deg = EQUIVOCAL_WITHIN_DEG,
    equivocal_between_deg = EQUIVOCAL_BETWEEN_DEG,
    equivocal_within_threshold = EQUIVOCAL_WITHIN_THRESHOLD,
    equivocal_between_threshold = EQUIVOCAL_BETWEEN_THRESHOLD,
    classification_rule = "Between if E > tan(52.5 deg); Within if E < tan(37.5 deg); otherwise Equivocal",
    mcse_definition = "MCSE = SD / sqrt(R)",
    k2_level5_rule = "r_R_het saved as sign-only / structurally degenerate for K = 2",
    stringsAsFactors = FALSE
  )
}

add_metadata_columns <- function(df) {
  df$equivocal_within_deg <- EQUIVOCAL_WITHIN_DEG
  df$equivocal_between_deg <- EQUIVOCAL_BETWEEN_DEG
  df$equivocal_within_threshold <- EQUIVOCAL_WITHIN_THRESHOLD
  df$equivocal_between_threshold <- EQUIVOCAL_BETWEEN_THRESHOLD
  df
}

summarize_metrics <- function(df, group_cols, metric_cols) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarize(
      n_rep = n(),
      across(all_of(metric_cols),
             list(mean = mean_or_na, sd = sd_or_na, mcse = mcse_or_na),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
}

# ============================================================================
# 1. Data Generation Functions
# ============================================================================

allocate_group_sizes <- function(n, K) {
  n_per <- rep(floor(n / K), K)
  n_per[K] <- n - sum(n_per[-K])
  n_per
}

generate_heteroscedastic <- function(n, p, K, mean_shift = 1.0,
                                     spread_corr = 0.8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_per <- allocate_group_sizes(n, K)

  group_means <- matrix(0, K, p)
  for (j in seq_len(p)) {
    group_means[, j] <- seq(-mean_shift, mean_shift, length.out = K) * (1 + 0.2 * (j - 1))
  }

  base_sd <- seq(0.5, 2.0, length.out = K)
  Sigma_noise <- spread_corr * matrix(1, p, p) + (1 - spread_corr) * diag(p)
  group_sds <- matrix(0, K, p)
  for (k in seq_len(K)) {
    noise <- mvrnorm(1, mu = rep(0, p), Sigma = Sigma_noise)
    group_sds[k, ] <- pmax(base_sd[k] + 0.3 * noise, 0.1)
  }

  X <- matrix(0, n, p)
  Y <- integer(n)
  idx <- 1
  for (k in seq_len(K)) {
    for (j in seq_len(p)) {
      X[idx:(idx + n_per[k] - 1), j] <- rnorm(n_per[k], group_means[k, j], group_sds[k, j])
    }
    Y[idx:(idx + n_per[k] - 1)] <- k
    idx <- idx + n_per[k]
  }
  list(X = X, Y = as.factor(Y))
}

generate_spread_only <- function(n, p, K, spread_range = c(0.5, 3.0),
                                 spread_corr = 0.9, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_per <- allocate_group_sizes(n, K)

  base_sd <- seq(spread_range[1], spread_range[2], length.out = K)
  Sigma_noise <- spread_corr * matrix(1, p, p) + (1 - spread_corr) * diag(p)
  group_sds <- matrix(0, K, p)
  for (k in seq_len(K)) {
    noise <- mvrnorm(1, mu = rep(0, p), Sigma = Sigma_noise)
    group_sds[k, ] <- pmax(base_sd[k] + 0.2 * noise, 0.1)
  }

  X <- matrix(0, n, p)
  Y <- integer(n)
  idx <- 1
  for (k in seq_len(K)) {
    for (j in seq_len(p)) {
      X[idx:(idx + n_per[k] - 1), j] <- rnorm(n_per[k], 0, group_sds[k, j])
    }
    Y[idx:(idx + n_per[k] - 1)] <- k
    idx <- idx + n_per[k]
  }
  list(X = X, Y = as.factor(Y))
}

generate_homoscedastic <- function(n, p, K, mean_shift = 1.0,
                                   common_sd = 1.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_per <- allocate_group_sizes(n, K)

  group_means <- matrix(0, K, p)
  for (j in seq_len(p)) {
    group_means[, j] <- seq(-mean_shift, mean_shift, length.out = K)
  }

  X <- matrix(0, n, p)
  Y <- integer(n)
  idx <- 1
  for (k in seq_len(K)) {
    X[idx:(idx + n_per[k] - 1), ] <- mvrnorm(
      n_per[k],
      group_means[k, ],
      common_sd^2 * diag(p)
    )
    Y[idx:(idx + n_per[k] - 1)] <- k
    idx <- idx + n_per[k]
  }
  list(X = X, Y = as.factor(Y))
}

generate_nonlinear_variance <- function(n, p, K, mean_shift = 2.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_per <- allocate_group_sizes(n, K)

  group_means <- matrix(0, K, p)
  group_sds <- matrix(0, K, p)
  for (j in seq_len(p)) {
    group_means[, j] <- seq(1, 1 + mean_shift, length.out = K)
    group_sds[, j] <- 0.5 * group_means[, j]
  }

  X <- matrix(0, n, p)
  Y <- integer(n)
  idx <- 1
  for (k in seq_len(K)) {
    for (j in seq_len(p)) {
      X[idx:(idx + n_per[k] - 1), j] <- rnorm(n_per[k], group_means[k, j], group_sds[k, j])
    }
    Y[idx:(idx + n_per[k] - 1)] <- k
    idx <- idx + n_per[k]
  }
  list(X = X, Y = as.factor(Y))
}

generate_complete_null <- function(n, p, K, common_sd = 1.0,
                                   common_corr = 0.35, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_per <- allocate_group_sizes(n, K)
  Sigma_common <- common_sd^2 *
    (common_corr * matrix(1, p, p) + (1 - common_corr) * diag(p))

  X <- matrix(0, n, p)
  Y <- integer(n)
  idx <- 1
  for (k in seq_len(K)) {
    X[idx:(idx + n_per[k] - 1), ] <- mvrnorm(n_per[k], rep(0, p), Sigma_common)
    Y[idx:(idx + n_per[k] - 1)] <- k
    idx <- idx + n_per[k]
  }
  list(X = X, Y = as.factor(Y))
}

# ============================================================================
# 2. Single-Replication Analysis
# ============================================================================

run_one_analysis <- function(X, Y) {
  result <- iwaba_full(X, Y)
  l1 <- result$level1
  l2 <- result$level2
  l3 <- result$level3
  l4 <- result$level4
  l5 <- result$level5
  eta <- l1$eta_results

  class_waba_explicit <- classify_e_explicit(eta$E_waba)
  class_int_explicit <- classify_e_explicit(eta$E_iwaba)

  level5_mode <- if (l5$is_k2_degenerate) {
    "k2_sign_only"
  } else if (l5$n_r_R_het_defined == 0) {
    "undefined"
  } else {
    "graded"
  }

  data.frame(
    info_gain = l1$info_gain,
    delta_tr = sum(l1$V_R),
    radius_prop = l1$radius_prop,
    mean_eta_B_waba = mean(eta$eta_B_waba),
    mean_eta_B_int = mean(eta$eta_B_iwaba),
    mean_delta_eta = mean(eta$delta_eta),
    mean_E_waba = mean(eta$E_waba),
    mean_E_int = mean(eta$E_iwaba),
    n_reclassified = sum(class_waba_explicit != class_int_explicit),
    n_between_waba = sum(class_waba_explicit == "Between"),
    n_between_int = sum(class_int_explicit == "Between"),
    n_equivocal_waba = sum(class_waba_explicit == "Equivocal"),
    n_equivocal_int = sum(class_int_explicit == "Equivocal"),
    n_within_waba = sum(class_waba_explicit == "Within"),
    n_within_int = sum(class_int_explicit == "Within"),
    mean_V_C = mean(l1$V_C),
    mean_V_R = mean(l1$V_R),
    mean_V_W = mean(l1$V_W),
    mean_abs_bw_int = l2$mean_abs_bw_int,
    mean_abs_ww_int = l2$mean_abs_ww_int,
    mean_abs_bw_waba = l2$mean_abs_bw_waba,
    mean_abs_ww_waba = l2$mean_abs_ww_waba,
    consistency_rate_int = l3$consistency_rate_int,
    consistency_rate_waba = l3$consistency_rate_waba,
    mean_pi_R = mean(l4$pi_R),
    mean_E_C_int = mean(l4$E_C_int),
    mean_E_R_int = mean(l4$E_R_int),
    mean_E_R_het_int = mean(l4$E_R_het_int),
    mean_V_R_het = mean(l4$V_R_het),
    mean_r_R_billard = l5$mean_r_R_billard,
    mean_r_R_het_signed = l5$mean_r_R_het_signed,
    prop_r_R_het_positive = l5$prop_r_R_het_positive,
    prop_r_R_het_negative = l5$prop_r_R_het_negative,
    prop_r_R_het_undefined = l5$prop_r_R_het_undefined,
    n_r_R_het_defined = l5$n_r_R_het_defined,
    n_r_R_het_undefined = l5$n_r_R_het_undefined,
    level5_mode = level5_mode,
    level5_k2_degenerate = l5$is_k2_degenerate,
    equivocal_within_deg = EQUIVOCAL_WITHIN_DEG,
    equivocal_between_deg = EQUIVOCAL_BETWEEN_DEG,
    equivocal_within_threshold = EQUIVOCAL_WITHIN_THRESHOLD,
    equivocal_between_threshold = EQUIVOCAL_BETWEEN_THRESHOLD,
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# 3. Main Simulation Runner
# ============================================================================

run_simulation <- function(n_reps = SIM_N_REPS,
                           K_values = SIM_K_VALUES,
                           n_values = SIM_N_VALUES,
                           p_values = SIM_P_VALUES,
                           scenarios = SIM_SCENARIOS,
                           verbose = TRUE) {
  p_values <- sort(unique(p_values))
  p_max <- max(p_values)

  generators <- list(
    heteroscedastic = generate_heteroscedastic,
    spread_only = generate_spread_only,
    homoscedastic = generate_homoscedastic,
    nonlinear_var = generate_nonlinear_variance,
    complete_null = generate_complete_null
  )

  total <- 0
  for (sc in scenarios) for (K in K_values) for (nn in n_values) {
    if (nn >= K * 5) for (pp in p_values) total <- total + 1
  }

  results <- vector("list", total * n_reps)
  ri <- 0
  counter <- 0

  for (scenario_index in seq_along(scenarios)) {
    scenario <- scenarios[scenario_index]
    gen_fun <- generators[[scenario]]

    for (K in K_values) {
      for (nn in n_values) {
        if (nn < K * 5) next

        for (rep in seq_len(n_reps)) {
          seed <- make_sim_seed(scenario_index, K, nn, rep)
          tryCatch({
            dat_template <- gen_fun(nn, p_max, K, seed = seed)

            for (pp in p_values) {
              if (rep == 1) {
                counter <- counter + 1
                if (verbose) {
                  cat(sprintf("[%d/%d] %s, K=%d, n=%d, p=%d (matched template)\n",
                              counter, total, scenario, K, nn, pp))
                }
              }

              dat <- subset_sim_data(dat_template, pp)
              metrics <- run_one_analysis(dat$X, dat$Y)
              metrics$scenario <- scenario
              metrics$K <- K
              metrics$n <- nn
              metrics$p <- pp
              metrics$rep <- rep
              ri <- ri + 1
              results[[ri]] <- metrics
            }
          }, error = function(e) {
            if (verbose) cat(sprintf("  Error in rep %d: %s\n", rep, e$message))
          })
        }
      }
    }
  }

  do.call(rbind, results[seq_len(ri)])
}

# ============================================================================
# 4. Comparator Experiment
# ============================================================================

welch_oneway_matrix <- function(X, g, varnames = NULL, alpha = 0.05) {
  X <- as.matrix(X)
  p <- ncol(X)
  if (is.null(varnames)) {
    varnames <- if (!is.null(colnames(X))) colnames(X) else paste0("V", seq_len(p))
  }

  rows <- vector("list", p)
  for (j in seq_len(p)) {
    p_value <- tryCatch(
      oneway.test(X[, j] ~ g, var.equal = FALSE)$p.value,
      error = function(e) NA_real_
    )
    rows[[j]] <- data.frame(
      variable = varnames[j],
      p_value = p_value,
      reject = !is.na(p_value) && p_value < alpha,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

run_one_comparator <- function(X, Y, alpha = COMPARATOR_ALPHA) {
  base_metrics <- run_one_analysis(X, Y)
  bf <- iwaba_brown_forsythe_matrix(X, Y, alpha = alpha)
  welch <- welch_oneway_matrix(X, Y, alpha = alpha)

  data.frame(
    mean_eta_B_waba = base_metrics$mean_eta_B_waba,
    mean_E_waba = base_metrics$mean_E_waba,
    n_between_waba = base_metrics$n_between_waba,
    any_between_waba = as.integer(base_metrics$n_between_waba > 0),
    n_equivocal_waba = base_metrics$n_equivocal_waba,
    mean_eta_B_int = base_metrics$mean_eta_B_int,
    mean_delta_eta = base_metrics$mean_delta_eta,
    delta_tr = base_metrics$delta_tr,
    mean_pi_R = base_metrics$mean_pi_R,
    mean_E_R_het_int = base_metrics$mean_E_R_het_int,
    info_gain = base_metrics$info_gain,
    bf_reject_count = sum(bf$reject, na.rm = TRUE),
    bf_reject_rate = mean(bf$reject, na.rm = TRUE),
    bf_any_reject = as.integer(sum(bf$reject, na.rm = TRUE) > 0),
    bf_min_p = min_or_na(bf$p_value),
    welch_reject_count = sum(welch$reject, na.rm = TRUE),
    welch_reject_rate = mean(welch$reject, na.rm = TRUE),
    welch_any_reject = as.integer(sum(welch$reject, na.rm = TRUE) > 0),
    welch_min_p = min_or_na(welch$p_value),
    stringsAsFactors = FALSE
  )
}

run_comparator_experiment <- function(n_reps = SIM_N_REPS,
                                      scenarios = COMPARATOR_SCENARIOS,
                                      K_values = SIM_K_VALUES,
                                      n = FOCAL_N,
                                      p = FOCAL_P,
                                      alpha = COMPARATOR_ALPHA,
                                      verbose = TRUE) {
  generators <- list(
    heteroscedastic = generate_heteroscedastic,
    spread_only = generate_spread_only,
    homoscedastic = generate_homoscedastic,
    nonlinear_var = generate_nonlinear_variance,
    complete_null = generate_complete_null
  )

  scenario_lookup <- setNames(seq_along(SIM_SCENARIOS), SIM_SCENARIOS)
  results <- vector("list", length(scenarios) * length(K_values) * n_reps)
  ri <- 0

  for (scenario in scenarios) {
    scenario_index <- scenario_lookup[[scenario]]
    gen_fun <- generators[[scenario]]

    for (K in K_values) {
      if (verbose) {
        cat(sprintf("Comparator: %s, K=%d, n=%d, p=%d\n", scenario, K, n, p))
      }
      for (rep in seq_len(n_reps)) {
        seed <- make_sim_seed(scenario_index, K, n, rep)
        dat <- gen_fun(n, p, K, seed = seed)
        metrics <- run_one_comparator(dat$X, dat$Y, alpha = alpha)
        metrics$scenario <- scenario
        metrics$K <- K
        metrics$n <- n
        metrics$p <- p
        metrics$rep <- rep
        ri <- ri + 1
        results[[ri]] <- metrics
      }
    }
  }

  raw_df <- do.call(rbind, results[seq_len(ri)])
  summary_df <- summarize_metrics(
    raw_df,
    group_cols = c("scenario", "K", "n", "p"),
    metric_cols = c(
      "mean_eta_B_waba", "mean_E_waba", "n_between_waba", "any_between_waba",
      "n_equivocal_waba", "mean_eta_B_int", "mean_delta_eta", "mean_pi_R",
      "mean_E_R_het_int", "info_gain", "delta_tr",
      "bf_reject_count", "bf_reject_rate", "bf_any_reject", "bf_min_p",
      "welch_reject_count", "welch_reject_rate", "welch_any_reject", "welch_min_p"
    )
  )
  summary_df$alpha <- alpha
  add_metadata_columns(summary_df)
}

# ============================================================================
# 5. Run Simulations
# ============================================================================

cat("================================================================\n")
cat("I-WABA Simulation Study v5\n")
cat("================================================================\n\n")

sim_results <- run_simulation(
  n_reps = SIM_N_REPS,
  K_values = SIM_K_VALUES,
  n_values = SIM_N_VALUES,
  p_values = SIM_P_VALUES,
  scenarios = SIM_SCENARIOS,
  verbose = TRUE
)

saveRDS(sim_results, file = RESULTS_RDS)
write.csv(sim_results, file = RESULTS_CSV, row.names = FALSE)

metadata_df <- simulation_metadata_df()
saveRDS(metadata_df, file = METADATA_RDS)
write.csv(metadata_df, file = METADATA_CSV, row.names = FALSE)

cat(sprintf("\nSaved %d raw rows to %s / %s\n", nrow(sim_results), RESULTS_RDS, RESULTS_CSV))

# ============================================================================
# 6. Summary Outputs
# ============================================================================

metric_cols_main <- c(
  "info_gain", "delta_tr", "radius_prop", "mean_eta_B_waba", "mean_eta_B_int",
  "mean_delta_eta", "n_reclassified", "mean_pi_R", "mean_E_R_het_int",
  "mean_r_R_billard", "mean_r_R_het_signed", "prop_r_R_het_positive",
  "prop_r_R_het_negative", "prop_r_R_het_undefined", "mean_abs_bw_int",
  "mean_abs_bw_waba", "consistency_rate_int", "consistency_rate_waba"
)

focal <- sim_results %>% filter(n == FOCAL_N, p == FOCAL_P)
summary_table <- summarize_metrics(focal, c("scenario", "K"), metric_cols_main) %>%
  left_join(
    focal %>%
      group_by(scenario, K) %>%
      summarize(
        level5_mode = summarize_level5_mode(level5_mode),
        n_r_R_het_defined_mean = mean_or_na(n_r_R_het_defined),
        n_r_R_het_undefined_mean = mean_or_na(n_r_R_het_undefined),
        .groups = "drop"
      ),
    by = c("scenario", "K")
  ) %>%
  rename(
    G_mean = info_gain_mean,
    G_sd = info_gain_sd,
    G_mcse = info_gain_mcse,
    Delta_tr = delta_tr_mean,
    Delta_tr_sd = delta_tr_sd,
    Delta_tr_mcse = delta_tr_mcse,
    radius_prop = radius_prop_mean,
    radius_prop_sd = radius_prop_sd,
    radius_prop_mcse = radius_prop_mcse,
    eta_B_waba = mean_eta_B_waba_mean,
    eta_B_waba_sd = mean_eta_B_waba_sd,
    eta_B_waba_mcse = mean_eta_B_waba_mcse,
    eta_B_int = mean_eta_B_int_mean,
    eta_B_int_sd = mean_eta_B_int_sd,
    eta_B_int_mcse = mean_eta_B_int_mcse,
    delta_eta = mean_delta_eta_mean,
    delta_eta_sd = mean_delta_eta_sd,
    delta_eta_mcse = mean_delta_eta_mcse,
    n_reclass = n_reclassified_mean,
    n_reclass_sd = n_reclassified_sd,
    n_reclass_mcse = n_reclassified_mcse,
    pi_R = mean_pi_R_mean,
    pi_R_sd = mean_pi_R_sd,
    pi_R_mcse = mean_pi_R_mcse,
    E_R_het = mean_E_R_het_int_mean,
    E_R_het_sd = mean_E_R_het_int_sd,
    E_R_het_mcse = mean_E_R_het_int_mcse,
    r_R_billard = mean_r_R_billard_mean,
    r_R_billard_sd = mean_r_R_billard_sd,
    r_R_billard_mcse = mean_r_R_billard_mcse,
    r_R_het_signed = mean_r_R_het_signed_mean,
    r_R_het_signed_sd = mean_r_R_het_signed_sd,
    r_R_het_signed_mcse = mean_r_R_het_signed_mcse,
    r_R_het_pos = prop_r_R_het_positive_mean,
    r_R_het_neg = prop_r_R_het_negative_mean,
    r_R_het_undef = prop_r_R_het_undefined_mean,
    bw_int = mean_abs_bw_int_mean,
    bw_int_sd = mean_abs_bw_int_sd,
    bw_int_mcse = mean_abs_bw_int_mcse,
    bw_waba = mean_abs_bw_waba_mean,
    bw_waba_sd = mean_abs_bw_waba_sd,
    bw_waba_mcse = mean_abs_bw_waba_mcse,
    cons_int = consistency_rate_int_mean,
    cons_int_sd = consistency_rate_int_sd,
    cons_int_mcse = consistency_rate_int_mcse,
    cons_waba = consistency_rate_waba_mean,
    cons_waba_sd = consistency_rate_waba_sd,
    cons_waba_mcse = consistency_rate_waba_mcse
  ) %>%
  mutate(
    level5_note = mapply(format_level5_note, level5_mode, r_R_het_pos, r_R_het_neg, r_R_het_undef),
    r_R_het_plot = if_else(level5_mode == "graded", r_R_het_signed, NA_real_)
  ) %>%
  add_metadata_columns()

saveRDS(summary_table, file = FOCAL_SUMMARY_RDS)
write.csv(summary_table, file = FOCAL_SUMMARY_CSV, row.names = FALSE)

p_sensitivity <- sim_results %>%
  filter(K == P_SENSITIVITY_K, n == P_SENSITIVITY_N) %>%
  summarize_metrics(c("scenario", "p"), c(
    "info_gain", "delta_tr", "radius_prop", "mean_eta_B_waba", "mean_eta_B_int",
    "mean_delta_eta", "mean_pi_R", "mean_E_R_het_int",
    "mean_r_R_billard", "mean_r_R_het_signed"
  )) %>%
  add_metadata_columns()

saveRDS(p_sensitivity, file = P_SENSITIVITY_RDS)
write.csv(p_sensitivity, file = P_SENSITIVITY_CSV, row.names = FALSE)

heatmap_data <- sim_results %>%
  filter(scenario == "heteroscedastic") %>%
  group_by(K, p) %>%
  summarize(
    n_rep = n(),
    G_mean = mean_or_na(info_gain),
    G_sd = sd_or_na(info_gain),
    G_mcse = mcse_or_na(info_gain),
    .groups = "drop"
  ) %>%
  add_metadata_columns()

saveRDS(heatmap_data, file = HEATMAP_RDS)
write.csv(heatmap_data, file = HEATMAP_CSV, row.names = FALSE)

complete_null_summary <- sim_results %>%
  filter(scenario == "complete_null") %>%
  summarize_metrics(c("K", "n", "p"), c(
    "info_gain", "delta_tr", "radius_prop", "mean_eta_B_waba", "mean_E_waba",
    "mean_eta_B_int", "mean_E_int", "mean_delta_eta", "mean_pi_R",
    "mean_E_R_het_int", "mean_r_R_billard", "mean_r_R_het_signed",
    "prop_r_R_het_positive", "prop_r_R_het_negative", "prop_r_R_het_undefined"
  )) %>%
  left_join(
    sim_results %>%
      filter(scenario == "complete_null") %>%
      group_by(K, n, p) %>%
      summarize(level5_mode = summarize_level5_mode(level5_mode), .groups = "drop"),
    by = c("K", "n", "p")
  ) %>%
  mutate(
    scenario = "complete_null",
    level5_note = mapply(format_level5_note, level5_mode,
                         prop_r_R_het_positive_mean,
                         prop_r_R_het_negative_mean,
                         prop_r_R_het_undefined_mean)
  ) %>%
  add_metadata_columns()

saveRDS(complete_null_summary, file = COMPLETE_NULL_RDS)
write.csv(complete_null_summary, file = COMPLETE_NULL_CSV, row.names = FALSE)

comparator_summary <- run_comparator_experiment()
saveRDS(comparator_summary, file = COMPARATOR_RDS)
write.csv(comparator_summary, file = COMPARATOR_CSV, row.names = FALSE)

# ============================================================================
# 7. Console Summaries
# ============================================================================

cat("\n================================================================\n")
cat("TABLE: Focal case n=500, p=10\n")
cat("================================================================\n")

cat("\n--- Level I: Entity Analysis ---\n")
cat(sprintf("%-18s %3s %8s %8s %9s %9s %6s %8s\n",
            "Scenario", "K", "G", "MCSE(G)", "RadProp",
            "eta_INT", "D.eta", "Reclass"))
for (i in seq_len(nrow(summary_table))) {
  r <- summary_table[i, ]
  cat(sprintf("%-18s %3d %8.3f %8.3f %9.3f %9.3f %6.3f %8.1f\n",
              r$scenario, r$K, r$G_mean, r$G_mcse, r$radius_prop,
              r$eta_B_int, r$delta_eta, r$n_reclass))
}

cat("\n--- Level IV: Dispersion-Source Diagnostics ---\n")
cat(sprintf("%-18s %3s %8s %10s %10s\n", "Scenario", "K", "pi_R", "E_R,het", "MCSE(E)"))
for (i in seq_len(nrow(summary_table))) {
  r <- summary_table[i, ]
  cat(sprintf("%-18s %3d %8.3f %10.3f %10.3f\n",
              r$scenario, r$K, r$pi_R, r$E_R_het, r$E_R_het_mcse))
}

cat("\n--- Level V: Dispersion Association ---\n")
cat(sprintf("%-18s %3s %10s %12s %s\n", "Scenario", "K", "r_R^Bill", "r_R,het", "Note"))
for (i in seq_len(nrow(summary_table))) {
  r <- summary_table[i, ]
  cat(sprintf("%-18s %3d %10s %12s %s\n",
              r$scenario, r$K,
              format_num_or_dash(r$r_R_billard),
              format_num_or_dash(r$r_R_het_signed),
              r$level5_note))
}

cat("\n--- Comparator Experiment (spread_only + complete_null; n=500, p=10) ---\n")
cat(sprintf("%-15s %3s %8s %8s %8s %8s %8s\n",
            "Scenario", "K", "WABA_B", "BF_any", "Welch", "I-eta", "E_Rhet"))
for (i in seq_len(nrow(comparator_summary))) {
  r <- comparator_summary[i, ]
  cat(sprintf("%-15s %3d %8.3f %8.3f %8.3f %8.3f %8.3f\n",
              r$scenario, r$K,
              r$any_between_waba_mean,
              r$bf_any_reject_mean,
              r$welch_any_reject_mean,
              r$mean_eta_B_int_mean,
              r$mean_E_R_het_int_mean))
}

# ============================================================================
# 8. Publication Figures
# ============================================================================

cat("\n================================================================\n")
cat("Generating publication figures...\n")
cat("================================================================\n")

scenario_colors <- c(
  "heteroscedastic" = "#1f77b4",
  "spread_only" = "#ff7f0e",
  "complete_null" = "#9467bd",
  "homoscedastic" = "#2ca02c",
  "nonlinear_var" = "#d62728"
)
scenario_labels <- c(
  "heteroscedastic" = "Heteroscedastic",
  "spread_only" = "Spread-Only",
  "complete_null" = "Complete Null",
  "homoscedastic" = "Homoscedastic",
  "nonlinear_var" = "Nonlinear-Var"
)

fig1_data <- summary_table %>%
  mutate(scenario_label = scenario_labels[scenario])

p1a <- ggplot(fig1_data %>% filter(scenario %in% c("heteroscedastic", "homoscedastic", "nonlinear_var")),
              aes(x = K, y = G_mean, color = scenario_label, shape = scenario_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(title = "(a) Mean-Structured Scenarios",
       x = "K", y = "Information Gain G", color = "Scenario", shape = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p1b <- ggplot(fig1_data %>% filter(scenario %in% c("spread_only", "complete_null")),
              aes(x = K, y = G_mean, color = scenario_label, shape = scenario_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(title = "(b) Spread-Only and Complete Null",
       x = "K", y = "Information Gain G", color = "Scenario", shape = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

fig1 <- grid.arrange(p1a, p1b, ncol = 2, widths = c(1.3, 1))
ggsave(FIG_INFO_GAIN_FILE, fig1, width = 12, height = 5)

eta_grid <- sim_results %>%
  filter(p == FOCAL_P, scenario %in% c("heteroscedastic", "spread_only", "complete_null")) %>%
  group_by(scenario, K, n) %>%
  summarize(delta_eta = mean_or_na(mean_delta_eta), .groups = "drop") %>%
  mutate(scenario_label = scenario_labels[scenario])

fig2 <- ggplot(eta_grid, aes(x = K, y = delta_eta, color = factor(n), shape = factor(n))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~scenario_label, scales = "free_y") +
  labs(title = "Improvement in Mean Between-Group Eta",
       x = "K", y = expression(bar(Delta) * eta[B]),
       color = "n", shape = "n") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(FIG_ETA_IMPROVEMENT_FILE, fig2, width = 12, height = 5)

p3a <- ggplot(fig1_data, aes(x = K, y = pi_R, color = scenario_label, shape = scenario_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(title = expression("(a) Mean Radius Share " * bar(pi)[R]),
       x = "K", y = expression(bar(pi)[R]),
       color = "Scenario", shape = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p3b <- ggplot(fig1_data, aes(x = K, y = E_R_het, color = scenario_label, shape = scenario_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(title = expression("(b) Heterogeneity Index " * bar(E)[R*","*het]^"(Int)"),
       x = "K", y = expression(bar(E)[R*","*het]^"(Int)"),
       color = "Scenario", shape = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

fig3 <- grid.arrange(p3a, p3b, ncol = 2)
ggsave(FIG_DISPERSION_FILE, fig3, width = 12, height = 5)

fig4 <- ggplot(heatmap_data, aes(x = factor(K), y = factor(p), fill = G_mean)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", G_mean)), size = 4) +
  scale_fill_gradient2(low = "white", mid = "steelblue", high = "darkblue",
                       midpoint = median(heatmap_data$G_mean)) +
  labs(title = "Information Gain (Heteroscedastic)",
       x = "K", y = "p", fill = "G") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")
ggsave(FIG_HEATMAP_FILE, fig4, width = 8, height = 6)

fig5 <- ggplot(fig1_data %>% filter(K > 2, !is.na(r_R_het_plot)),
               aes(x = K, y = r_R_het_plot, color = scenario_label, shape = scenario_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(title = expression("Level V: Signed Heterogeneity-Based Radius Correlation " * bar(r)[R*","*het]),
       subtitle = "K = 2 is sign-only; complete null omitted where centered radii are undefined",
       x = "K", y = expression(bar(r)[R*","*het]),
       color = "Scenario", shape = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(FIG_RADIUS_FILE, fig5, width = 8, height = 5)

# ============================================================================
# 9. LaTeX Tables
# ============================================================================

cat("\n================================================================\n")
cat("LaTeX Tables\n")
cat("================================================================\n\n")

cat("% Table 3: Level I Entity Analysis\n")
cat("\\begin{table}[ht]\n\\centering\n")
cat("\\caption{Simulation results: Level~I entity analysis ($n=500$, $p=10$, 100 replications). The equivocal band is defined explicitly by $\\tan(37.5^\\circ) < E < \\tan(52.5^\\circ)$.}\n")
cat("\\label{tab:sim_results}\n")
cat("\\resizebox{\\textwidth}{!}{%\n")
cat("\\begin{tabular}{llcccccc}\n\\toprule\n")
cat("Scenario & $K$ & $\\mathcal{G}$ & MCSE$(\\mathcal{G})$ & $\\bar{\\pi}_R$ & $\\bar{\\eta}_B^{\\text{WABA}}$ & $\\bar{\\eta}_B^{\\text{I-WABA}}$ & Reclass. \\\\\n")
cat("\\midrule\n")

prev_sc <- ""
for (i in seq_len(nrow(summary_table))) {
  r <- summary_table[i, ]
  sc <- scenario_labels[r$scenario]
  if (prev_sc != "" && sc != prev_sc) cat("\\midrule\n")
  cat(sprintf("%s & %d & %.3f & %.3f & %.3f & %.3f & %.3f & %.1f \\\\\n",
              sc, r$K, r$G_mean, r$G_mcse, r$radius_prop,
              r$eta_B_waba, r$eta_B_int, r$n_reclass))
  prev_sc <- sc
}
cat("\\bottomrule\n\\end{tabular}%\n}\n\\end{table}\n\n")

cat("% Level IV: Dispersion-Source Diagnostics\n")
cat("\\begin{table}[ht]\n\\centering\n")
cat("\\caption{Level~IV dispersion-source diagnostics ($n=500$, $p=10$). MCSEs are available in the saved summary output.}\n")
cat("\\label{tab:sim_level_iv}\n")
cat("\\begin{tabular}{llcc}\n\\toprule\n")
cat("Scenario & $K$ & $\\bar{\\pi}_R$ & $\\bar{E}_{R,\\text{het}}^{(\\text{Int})}$ \\\\\n")
cat("\\midrule\n")

prev_sc <- ""
for (i in seq_len(nrow(summary_table))) {
  r <- summary_table[i, ]
  sc <- scenario_labels[r$scenario]
  if (prev_sc != "" && sc != prev_sc) cat("\\midrule\n")
  cat(sprintf("%s & %d & %.3f & %.3f \\\\\n", sc, r$K, r$pi_R, r$E_R_het))
  prev_sc <- sc
}
cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")

cat("% Level V: Dispersion Association\n")
cat("\\begin{table}[ht]\n\\centering\n")
cat("\\caption{Level~V dispersion association ($n=500$, $p=10$). $\\bar{r}_R^{(\\text{Billard})}$ is reported on its natural $[0,1]$ scale. For $K=2$, $r_{R,\\text{het}}$ is saved as a sign-only / structurally degenerate diagnostic. When centered radii vanish, $r_{R,\\text{het}}$ is undefined.}\n")
cat("\\label{tab:sim_level_v}\n")
cat("\\begin{tabular}{llccc}\n\\toprule\n")
cat("Scenario & $K$ & $\\bar{r}_R^{(\\text{Billard})}$ & $\\bar{r}_{R,\\text{het}}$ & Note \\\\\n")
cat("\\midrule\n")

prev_sc <- ""
for (i in seq_len(nrow(summary_table))) {
  r <- summary_table[i, ]
  sc <- scenario_labels[r$scenario]
  if (prev_sc != "" && sc != prev_sc) cat("\\midrule\n")
  note_text <- if (nzchar(r$level5_note)) r$level5_note else "--"
  cat(sprintf("%s & %d & %s & %s & %s \\\\\n",
              sc, r$K,
              format_num_or_dash(r$r_R_billard),
              format_num_or_dash(r$r_R_het_signed),
              note_text))
  prev_sc <- sc
}
cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

cat("\n================================================================\n")
cat("Saved outputs:\n")
cat(sprintf("  %s\n", RESULTS_RDS))
cat(sprintf("  %s\n", RESULTS_CSV))
cat(sprintf("  %s\n", METADATA_RDS))
cat(sprintf("  %s\n", METADATA_CSV))
cat(sprintf("  %s\n", FOCAL_SUMMARY_CSV))
cat(sprintf("  %s\n", P_SENSITIVITY_CSV))
cat(sprintf("  %s\n", HEATMAP_CSV))
cat(sprintf("  %s\n", COMPLETE_NULL_CSV))
cat(sprintf("  %s\n", COMPARATOR_CSV))
cat(sprintf("  %s\n", FIG_INFO_GAIN_FILE))
cat(sprintf("  %s\n", FIG_ETA_IMPROVEMENT_FILE))
cat(sprintf("  %s\n", FIG_DISPERSION_FILE))
cat(sprintf("  %s\n", FIG_HEATMAP_FILE))
cat(sprintf("  %s\n", FIG_RADIUS_FILE))
cat("================================================================\n")
cat("Simulation study complete.\n")
cat("================================================================\n")
