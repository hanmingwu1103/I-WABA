################################################################################
## Correlation-ordering sensitivity summary for I-WABA
## Purpose: quantify how often r_Between^(Int) < r_Between in the focal
## simulation design used in the manuscript.
################################################################################

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required.")
}

source("iwaba_functions.R")

RESULTS_DIR <- "results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

OUTPUT_VERSION <- "v1"
DETAIL_CSV <- file.path(
  RESULTS_DIR,
  sprintf("simulation_correlation_ordering_results_%s.csv", OUTPUT_VERSION)
)
SUMMARY_CSV <- file.path(
  RESULTS_DIR,
  sprintf("simulation_correlation_ordering_summary_%s.csv", OUTPUT_VERSION)
)

SIM_N_REPS <- 100
SIM_K_VALUES <- c(2, 3, 5, 10)
SIM_SCENARIOS <- c(
  "heteroscedastic",
  "spread_only",
  "homoscedastic",
  "nonlinear_var",
  "complete_null"
)
FOCAL_N <- 500
FOCAL_P <- 10

allocate_group_sizes <- function(n, K) {
  n_per <- rep(floor(n / K), K)
  n_per[K] <- n - sum(n_per[-K])
  n_per
}

make_sim_seed <- function(scenario_index, K, n, rep) {
  scenario_index * 1000000 + K * 10000 + n * 10 + rep
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
    noise <- MASS::mvrnorm(1, mu = rep(0, p), Sigma = Sigma_noise)
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
    noise <- MASS::mvrnorm(1, mu = rep(0, p), Sigma = Sigma_noise)
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
    X[idx:(idx + n_per[k] - 1), ] <- MASS::mvrnorm(
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
    X[idx:(idx + n_per[k] - 1), ] <- MASS::mvrnorm(n_per[k], rep(0, p), Sigma_common)
    Y[idx:(idx + n_per[k] - 1)] <- k
    idx <- idx + n_per[k]
  }
  list(X = X, Y = as.factor(Y))
}

generators <- list(
  heteroscedastic = generate_heteroscedastic,
  spread_only = generate_spread_only,
  homoscedastic = generate_homoscedastic,
  nonlinear_var = generate_nonlinear_variance,
  complete_null = generate_complete_null
)

one_rep_summary <- function(X, Y) {
  fit <- iwaba_full(X, Y)
  l2 <- fit$level2
  p <- ncol(X)
  pairs <- combn(seq_len(p), 2)

  r_int <- vapply(
    seq_len(ncol(pairs)),
    function(idx) l2$r_between_int[pairs[1, idx], pairs[2, idx]],
    numeric(1)
  )
  r_classical <- vapply(
    seq_len(ncol(pairs)),
    function(idx) l2$r_between_classical[pairs[1, idx], pairs[2, idx]],
    numeric(1)
  )

  ok <- is.finite(r_int) & is.finite(r_classical)
  diffs <- r_int[ok] - r_classical[ok]
  n_pairs <- length(diffs)
  if (n_pairs == 0) {
    return(data.frame(
      n_pairs = 0,
      prop_corr_decrease = NA_real_,
      prop_corr_increase = NA_real_,
      prop_corr_equal = NA_real_,
      mean_signed_diff = NA_real_,
      mean_negative_diff = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  neg <- diffs < -1e-12
  pos <- diffs > 1e-12
  eq <- !(neg | pos)

  data.frame(
    n_pairs = n_pairs,
    prop_corr_decrease = mean(neg),
    prop_corr_increase = mean(pos),
    prop_corr_equal = mean(eq),
    mean_signed_diff = mean(diffs),
    mean_negative_diff = if (any(neg)) mean(diffs[neg]) else 0,
    stringsAsFactors = FALSE
  )
}

results <- vector("list", length(SIM_SCENARIOS) * length(SIM_K_VALUES) * SIM_N_REPS)
idx <- 0

for (scenario_index in seq_along(SIM_SCENARIOS)) {
  scenario <- SIM_SCENARIOS[scenario_index]
  generator <- generators[[scenario]]
  for (K in SIM_K_VALUES) {
    for (rep in seq_len(SIM_N_REPS)) {
      seed <- make_sim_seed(scenario_index, K, FOCAL_N, rep)
      dat <- generator(n = FOCAL_N, p = FOCAL_P, K = K, seed = seed)
      smry <- one_rep_summary(dat$X, dat$Y)
      idx <- idx + 1
      results[[idx]] <- data.frame(
        scenario = scenario,
        K = K,
        n = FOCAL_N,
        p = FOCAL_P,
        rep = rep,
        smry,
        stringsAsFactors = FALSE
      )
    }
  }
}

detail_df <- do.call(rbind, results)
summary_df <- aggregate(
  cbind(
    prop_corr_decrease,
    prop_corr_increase,
    prop_corr_equal,
    mean_signed_diff,
    mean_negative_diff
  ) ~ scenario + K + n + p,
  data = detail_df,
  FUN = mean
)

summary_df$mc_reps <- SIM_N_REPS
summary_df$total_pairs <- FOCAL_P * (FOCAL_P - 1) / 2
summary_df <- summary_df[, c(
  "scenario", "K", "n", "p", "mc_reps", "total_pairs",
  "prop_corr_decrease", "prop_corr_increase", "prop_corr_equal",
  "mean_signed_diff", "mean_negative_diff"
)]

write.csv(detail_df, DETAIL_CSV, row.names = FALSE)
write.csv(summary_df, SUMMARY_CSV, row.names = FALSE)

print(summary_df)
