################################################################################
## I-WABA Inferential Simulation Study (v1)
## Title: Interval-based Within-And-Between Analysis
## Author: Han-Ming Wu (National Chengchi University)
##
## Revised Section 3.6 design:
##   - K = 5, p = 5
##   - n in {100, 500}, equal group sizes
##   - Monte Carlo replications R = 300
##   - Within-group bootstrap resamples B = 399
##
## Targets:
##   - Coverage and mean interval length for:
##       delta_1^(E)  = V_C^(1) + V_R^(1) - V_W^(1)
##       delta_12^(A) = (D_B,12^(Int))^2 - (D_W,12^(Int))^2
##       r_R,het(X1, X2)
##   - Size / power for:
##       Brown-Forsythe equal-variance test on X1
##       Percentile-bootstrap CI exclusion test for H0: r_R,het(X1, X2) = 0
##
## Output:
##   - simulation_inference_results_v1.rds / .csv
##   - simulation_inference_coverage_table_v1.rds / .csv
##   - simulation_inference_size_power_table_v1.rds / .csv
################################################################################

source("iwaba_functions.R")

OUTPUT_VERSION <- "v1"
RESULTS_RDS <- sprintf("simulation_inference_results_%s.rds", OUTPUT_VERSION)
RESULTS_CSV <- sprintf("simulation_inference_results_%s.csv", OUTPUT_VERSION)
COVERAGE_RDS <- sprintf("simulation_inference_coverage_table_%s.rds", OUTPUT_VERSION)
COVERAGE_CSV <- sprintf("simulation_inference_coverage_table_%s.csv", OUTPUT_VERSION)
SIZE_POWER_RDS <- sprintf("simulation_inference_size_power_table_%s.rds", OUTPUT_VERSION)
SIZE_POWER_CSV <- sprintf("simulation_inference_size_power_table_%s.csv", OUTPUT_VERSION)

DEFAULT_R <- 300L
DEFAULT_B <- 399L
ALPHA <- 0.05
K_FIXED <- 5L
P_FIXED <- 5L
N_VALUES <- c(100L, 500L)
SCENARIOS <- c("aligned_heteroscedastic", "homoscedastic_null", "orthogonal_heterogeneity")

parse_arg_value <- function(args, prefix, default) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) == 0) {
    default
  } else {
    sub(paste0("^", prefix), "", hit[1])
  }
}

args <- commandArgs(trailingOnly = TRUE)
MC_REPS <- as.integer(parse_arg_value(args, "--R=", DEFAULT_R))
BOOT_REPS <- as.integer(parse_arg_value(args, "--B=", DEFAULT_B))

if (is.na(MC_REPS) || MC_REPS <= 0L) stop("MC replications must be a positive integer.")
if (is.na(BOOT_REPS) || BOOT_REPS <= 0L) stop("Bootstrap resamples must be a positive integer.")

make_seed <- function(scenario_index, n, rep) {
  scenario_index * 1000000L + n * 1000L + rep
}

mean_pattern <- c(-1.00, -0.50, 0.00, 0.50, 1.00)
mean_scales <- c(1.40, 1.60, 1.45, 1.80, 1.35)

make_scenario_template <- function(scenario) {
  means <- outer(mean_pattern, mean_scales)
  colnames(means) <- paste0("X", seq_len(P_FIXED))

  if (scenario == "aligned_heteroscedastic") {
    pattern_x1 <- c(-0.36, -0.18, 0.00, 0.18, 0.36)
    pattern_x2 <- c(-0.24, -0.18, 0.12, 0.12, 0.18)
    sds <- cbind(
      0.95 + pattern_x1,
      1.05 + pattern_x2,
      1.10 + 0.80 * pattern_x1,
      1.00 + 0.60 * pattern_x2,
      1.15 + 0.70 * pattern_x1
    )
  } else if (scenario == "homoscedastic_null") {
    sd_base <- c(0.95, 1.05, 1.10, 1.00, 1.15)
    sds <- matrix(sd_base, nrow = K_FIXED, ncol = P_FIXED, byrow = TRUE)
  } else if (scenario == "orthogonal_heterogeneity") {
    pattern_x1 <- c(-0.36, -0.18, 0.00, 0.18, 0.36)
    pattern_x2 <- c(0.18, -0.36, 0.00, 0.36, -0.18)
    sds <- cbind(
      1.00 + pattern_x1,
      1.05 + pattern_x2,
      0.95 + 0.90 * pattern_x1,
      1.10 - 0.70 * pattern_x1,
      1.00 + 0.80 * pattern_x2
    )
  } else {
    stop("Unknown scenario: ", scenario)
  }

  if (any(sds <= 0)) {
    stop("Scenario ", scenario, " generated non-positive SDs.")
  }

  colnames(sds) <- paste0("X", seq_len(P_FIXED))
  rownames(means) <- paste0("G", seq_len(K_FIXED))
  rownames(sds) <- paste0("G", seq_len(K_FIXED))

  list(
    scenario = scenario,
    means = means,
    sds = sds,
    w_k = rep(1 / K_FIXED, K_FIXED)
  )
}

cov_to_corr_safe <- function(S) {
  if (all(dim(S) == c(1L, 1L))) {
    return(matrix(1, 1, 1))
  }
  cov_to_corr(S)
}

weighted_center <- function(M, weights) {
  colMeans_weighted <- colSums(sweep(M, 1, weights, "*"))
  sweep(M, 2, colMeans_weighted)
}

weighted_corr_pair <- function(x, y, weights, eps = 1e-15) {
  x_centered <- x - sum(weights * x)
  y_centered <- y - sum(weights * y)
  denom <- sqrt(sum(weights * x_centered^2) * sum(weights * y_centered^2))
  if (denom <= eps) {
    NA_real_
  } else {
    max(min(sum(weights * x_centered * y_centered) / denom, 1), -1)
  }
}

true_population_targets <- function(template) {
  weights <- template$w_k
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
  Cov_aug <- Cov_between_int + Cov_within

  eta_B_int <- sqrt(safe_div(diag(Cov_between_int), diag(Cov_aug)))
  eta_W_int <- sqrt(safe_div(diag(Cov_within), diag(Cov_aug)))
  r_between_int <- cov_to_corr_safe(Cov_between_int)
  r_within <- cov_to_corr_safe(Cov_within)

  D_B_12 <- eta_B_int[1] * eta_B_int[2] * r_between_int[1, 2]
  D_W_12 <- eta_W_int[1] * eta_W_int[2] * r_within[1, 2]

  list(
    delta_E_1 = V_C[1] + V_R[1] - V_W[1],
    delta_A_12 = D_B_12^2 - D_W_12^2,
    r_R_het_12 = true_rhet,
    r_R_het_defined = !is.na(true_rhet),
    bf_equal_var_null = all(abs(sds[, 1] - sds[1, 1]) < 1e-12),
    rhet_zero_null = !is.na(true_rhet) && abs(true_rhet) < 1e-12
  )
}

simulate_from_template <- function(template, n, seed) {
  if (n %% K_FIXED != 0L) {
    stop("n must be divisible by K for equal group sizes.")
  }
  set.seed(seed)

  n_per_group <- n %/% K_FIXED
  X <- matrix(NA_real_, nrow = n, ncol = P_FIXED)
  Y <- factor(rep(seq_len(K_FIXED), each = n_per_group), levels = seq_len(K_FIXED))
  colnames(X) <- paste0("X", seq_len(P_FIXED))

  row_start <- 1L
  for (k in seq_len(K_FIXED)) {
    row_end <- row_start + n_per_group - 1L
    Z <- matrix(rnorm(n_per_group * P_FIXED), nrow = n_per_group, ncol = P_FIXED)
    X[row_start:row_end, ] <- sweep(Z, 2, template$sds[k, ], "*")
    X[row_start:row_end, ] <- sweep(X[row_start:row_end, , drop = FALSE], 2, template$means[k, ], "+")
    row_start <- row_end + 1L
  }

  list(X = X, Y = Y)
}

find_pair_index <- function(pairs, i, j) {
  which(pairs[1, ] == i & pairs[2, ] == j)
}

extract_targets <- function(analysis) {
  pair_idx <- find_pair_index(analysis$level2$pairs, 1, 2)
  if (length(pair_idx) != 1L) {
    stop("Could not locate pair (1, 2) in Level II output.")
  }

  list(
    delta_E_1 = analysis$level1$V_C[1] + analysis$level1$V_R[1] - analysis$level1$V_W[1],
    delta_A_12 = analysis$level2$bw_int[pair_idx]^2 - analysis$level2$ww_int[pair_idx]^2,
    r_R_het_12 = analysis$level5$r_R_het[pair_idx],
    r_R_het_defined = analysis$level5$r_R_het_defined[pair_idx]
  )
}

brown_forsythe_pvalue <- function(x, g) {
  g <- factor(g)
  K <- nlevels(g)
  n <- length(x)
  group_medians <- tapply(x, g, median)
  deviations <- abs(x - group_medians[g])
  group_sizes <- as.numeric(table(g))
  group_means <- as.numeric(tapply(deviations, g, mean))
  grand_mean <- mean(deviations)

  ss_between <- sum(group_sizes * (group_means - grand_mean)^2)
  ss_within <- sum(tapply(deviations, g, function(z) sum((z - mean(z))^2)))

  df1 <- K - 1L
  df2 <- n - K
  if (df1 <= 0L || df2 <= 0L || ss_within <= 0) {
    return(NA_real_)
  }

  f_stat <- (ss_between / df1) / (ss_within / df2)
  1 - pf(f_stat, df1, df2)
}

percentile_interval <- function(x, alpha = ALPHA) {
  x <- x[is.finite(x)]
  if (length(x) < 25L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  probs <- c(alpha / 2, 1 - alpha / 2)
  stats::quantile(x, probs = probs, na.rm = TRUE, names = TRUE, type = 7)
}

bootstrap_targets <- function(X, Y, B) {
  group_indices <- split(seq_len(nrow(X)), Y)
  boot_delta_E <- numeric(B)
  boot_delta_A <- numeric(B)
  boot_rhet <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    boot_idx <- unlist(lapply(group_indices, function(idx) sample(idx, length(idx), replace = TRUE)),
      use.names = FALSE
    )
    analysis_b <- iwaba_full(X[boot_idx, , drop = FALSE], Y[boot_idx])
    targets_b <- extract_targets(analysis_b)
    boot_delta_E[b] <- targets_b$delta_E_1
    boot_delta_A[b] <- targets_b$delta_A_12
    if (isTRUE(targets_b$r_R_het_defined)) {
      boot_rhet[b] <- targets_b$r_R_het_12
    }
  }

  list(
    delta_E = boot_delta_E,
    delta_A = boot_delta_A,
    r_R_het = boot_rhet
  )
}

coverage_indicator <- function(ci, truth) {
  if (is.na(truth) || anyNA(ci)) {
    NA
  } else {
    truth >= ci[1] && truth <= ci[2]
  }
}

interval_length <- function(ci) {
  if (anyNA(ci)) NA_real_ else unname(ci[2] - ci[1])
}

ci_exclusion_reject <- function(ci, null_value = 0) {
  if (anyNA(ci)) {
    NA
  } else {
    (null_value < ci[1]) || (null_value > ci[2])
  }
}

format_cov_entry <- function(coverage, length_value) {
  if (is.na(coverage) || is.na(length_value)) {
    "--"
  } else {
    sprintf("%.3f [%.3f]", coverage, length_value)
  }
}

format_test_entry <- function(rate, role) {
  if (is.na(rate)) {
    "--"
  } else {
    sprintf("%.3f (%s)", rate, role)
  }
}

run_inference_simulation <- function(mc_reps, boot_reps) {
  templates <- lapply(SCENARIOS, make_scenario_template)
  true_targets <- lapply(templates, true_population_targets)
  names(templates) <- SCENARIOS
  names(true_targets) <- SCENARIOS

  n_total <- length(SCENARIOS) * length(N_VALUES) * mc_reps
  pb <- txtProgressBar(min = 0, max = n_total, style = 3)
  progress <- 0L

  results <- vector("list", n_total)
  row_id <- 1L

  for (scenario_index in seq_along(SCENARIOS)) {
    scenario_name <- SCENARIOS[scenario_index]
    template <- templates[[scenario_name]]
    truth <- true_targets[[scenario_name]]

    for (n in N_VALUES) {
      for (rep in seq_len(mc_reps)) {
        data_seed <- make_seed(scenario_index, n, rep)
        sim_data <- simulate_from_template(template, n = n, seed = data_seed)
        analysis <- iwaba_full(sim_data$X, sim_data$Y)
        targets <- extract_targets(analysis)
        boot_stats <- bootstrap_targets(sim_data$X, sim_data$Y, B = boot_reps)

        ci_delta_E <- percentile_interval(boot_stats$delta_E)
        ci_delta_A <- percentile_interval(boot_stats$delta_A)
        ci_rhet <- percentile_interval(boot_stats$r_R_het)

        bf_pvalue <- brown_forsythe_pvalue(sim_data$X[, 1], sim_data$Y)

        results[[row_id]] <- data.frame(
          scenario = scenario_name,
          n = n,
          rep = rep,
          true_delta_E_1 = truth$delta_E_1,
          true_delta_A_12 = truth$delta_A_12,
          true_r_R_het_12 = truth$r_R_het_12,
          true_r_R_het_defined = truth$r_R_het_defined,
          sample_delta_E_1 = targets$delta_E_1,
          sample_delta_A_12 = targets$delta_A_12,
          sample_r_R_het_12 = targets$r_R_het_12,
          sample_r_R_het_defined = targets$r_R_het_defined,
          ci_delta_E_lower = unname(ci_delta_E[1]),
          ci_delta_E_upper = unname(ci_delta_E[2]),
          ci_delta_A_lower = unname(ci_delta_A[1]),
          ci_delta_A_upper = unname(ci_delta_A[2]),
          ci_r_R_het_lower = unname(ci_rhet[1]),
          ci_r_R_het_upper = unname(ci_rhet[2]),
          cover_delta_E_1 = coverage_indicator(ci_delta_E, truth$delta_E_1),
          cover_delta_A_12 = coverage_indicator(ci_delta_A, truth$delta_A_12),
          cover_r_R_het_12 = coverage_indicator(ci_rhet, truth$r_R_het_12),
          len_delta_E_1 = interval_length(ci_delta_E),
          len_delta_A_12 = interval_length(ci_delta_A),
          len_r_R_het_12 = interval_length(ci_rhet),
          bf_pvalue_x1 = bf_pvalue,
          reject_bf_x1 = !is.na(bf_pvalue) && bf_pvalue < ALPHA,
          reject_r_R_het_zero = ci_exclusion_reject(ci_rhet, null_value = 0),
          boot_r_R_het_defined_rate = mean(is.finite(boot_stats$r_R_het)),
          stringsAsFactors = FALSE
        )

        row_id <- row_id + 1L
        progress <- progress + 1L
        setTxtProgressBar(pb, progress)
      }
    }
  }

  close(pb)

  do.call(rbind, results)
}

build_summary_tables <- function(results) {
  scenarios_in_order <- SCENARIOS
  summary_grid <- expand.grid(
    scenario = scenarios_in_order,
    n = N_VALUES,
    stringsAsFactors = FALSE
  )

  coverage_rows <- vector("list", nrow(summary_grid))
  size_rows <- vector("list", nrow(summary_grid))

  for (i in seq_len(nrow(summary_grid))) {
    scenario_name <- summary_grid$scenario[i]
    n_value <- summary_grid$n[i]
    subset_rows <- results[results$scenario == scenario_name & results$n == n_value, , drop = FALSE]

    true_rhet_defined <- unique(subset_rows$true_r_R_het_defined)
    if (length(true_rhet_defined) != 1L) stop("Unexpected heterogeneity target definition mismatch.")
    true_rhet_defined <- isTRUE(true_rhet_defined)

    bf_role <- if (scenario_name == "homoscedastic_null") "size" else "power"
    rhet_role <- if (!true_rhet_defined) {
      "undefined"
    } else if (abs(unique(subset_rows$true_r_R_het_12)) < 1e-12) {
      "size"
    } else {
      "power"
    }

    coverage_rows[[i]] <- data.frame(
      scenario = scenario_name,
      n = n_value,
      R = nrow(subset_rows),
      B = BOOT_REPS,
      true_delta_E_1 = unique(subset_rows$true_delta_E_1),
      coverage_delta_E_1 = mean(subset_rows$cover_delta_E_1, na.rm = TRUE),
      mean_ci_len_delta_E_1 = mean(subset_rows$len_delta_E_1, na.rm = TRUE),
      true_delta_A_12 = unique(subset_rows$true_delta_A_12),
      coverage_delta_A_12 = mean(subset_rows$cover_delta_A_12, na.rm = TRUE),
      mean_ci_len_delta_A_12 = mean(subset_rows$len_delta_A_12, na.rm = TRUE),
      true_r_R_het_12 = unique(subset_rows$true_r_R_het_12),
      true_r_R_het_defined = true_rhet_defined,
      coverage_r_R_het_12 = if (true_rhet_defined) mean(subset_rows$cover_r_R_het_12, na.rm = TRUE) else NA_real_,
      mean_ci_len_r_R_het_12 = if (true_rhet_defined) mean(subset_rows$len_r_R_het_12, na.rm = TRUE) else NA_real_,
      mean_boot_r_R_het_defined_rate = mean(subset_rows$boot_r_R_het_defined_rate, na.rm = TRUE),
      delta_E_display = format_cov_entry(
        mean(subset_rows$cover_delta_E_1, na.rm = TRUE),
        mean(subset_rows$len_delta_E_1, na.rm = TRUE)
      ),
      delta_A_display = format_cov_entry(
        mean(subset_rows$cover_delta_A_12, na.rm = TRUE),
        mean(subset_rows$len_delta_A_12, na.rm = TRUE)
      ),
      r_R_het_display = if (true_rhet_defined) {
        format_cov_entry(
          mean(subset_rows$cover_r_R_het_12, na.rm = TRUE),
          mean(subset_rows$len_r_R_het_12, na.rm = TRUE)
        )
      } else {
        "undefined at population level"
      },
      stringsAsFactors = FALSE
    )

    size_rows[[i]] <- data.frame(
      scenario = scenario_name,
      n = n_value,
      R = nrow(subset_rows),
      B = BOOT_REPS,
      bf_role = bf_role,
      bf_reject_rate_x1 = mean(subset_rows$reject_bf_x1, na.rm = TRUE),
      true_r_R_het_12 = unique(subset_rows$true_r_R_het_12),
      true_r_R_het_defined = true_rhet_defined,
      r_R_het_role = rhet_role,
      r_R_het_reject_rate = if (true_rhet_defined) mean(subset_rows$reject_r_R_het_zero, na.rm = TRUE) else NA_real_,
      mean_boot_r_R_het_defined_rate = mean(subset_rows$boot_r_R_het_defined_rate, na.rm = TRUE),
      bf_display = format_test_entry(mean(subset_rows$reject_bf_x1, na.rm = TRUE), bf_role),
      r_R_het_display = if (true_rhet_defined) {
        format_test_entry(mean(subset_rows$reject_r_R_het_zero, na.rm = TRUE), rhet_role)
      } else {
        "undefined at population level"
      },
      stringsAsFactors = FALSE
    )
  }

  list(
    coverage = do.call(rbind, coverage_rows),
    size_power = do.call(rbind, size_rows)
  )
}

cat("Running inferential simulation with:\n")
cat(sprintf("  K = %d, p = %d, n in {%s}\n", K_FIXED, P_FIXED, paste(N_VALUES, collapse = ", ")))
cat(sprintf("  Monte Carlo replications R = %d\n", MC_REPS))
cat(sprintf("  Bootstrap resamples B = %d\n", BOOT_REPS))
cat(sprintf("  Output version = %s\n\n", OUTPUT_VERSION))

results <- run_inference_simulation(MC_REPS, BOOT_REPS)
tables <- build_summary_tables(results)

saveRDS(results, RESULTS_RDS)
write.csv(results, RESULTS_CSV, row.names = FALSE)
saveRDS(tables$coverage, COVERAGE_RDS)
write.csv(tables$coverage, COVERAGE_CSV, row.names = FALSE)
saveRDS(tables$size_power, SIZE_POWER_RDS)
write.csv(tables$size_power, SIZE_POWER_CSV, row.names = FALSE)

cat("\nCoverage table:\n")
print(tables$coverage[, c("scenario", "n", "delta_E_display", "delta_A_display", "r_R_het_display")], row.names = FALSE)

cat("\nSize / power table:\n")
print(tables$size_power[, c("scenario", "n", "bf_display", "r_R_het_display")], row.names = FALSE)

cat("\nSaved outputs:\n")
cat(sprintf("  %s\n", RESULTS_RDS))
cat(sprintf("  %s\n", RESULTS_CSV))
cat(sprintf("  %s\n", COVERAGE_RDS))
cat(sprintf("  %s\n", COVERAGE_CSV))
cat(sprintf("  %s\n", SIZE_POWER_RDS))
cat(sprintf("  %s\n", SIZE_POWER_CSV))
