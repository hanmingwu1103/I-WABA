################################################################################
## I-WABA Simulation Study
## Title: Interval-based Within-And-Between Analysis
## Author: Han-Ming Wu
## Description: Simulation studies comparing Classical WABA and I-WABA
##              across varying K (number of classes), n, and p
################################################################################

# ============================================================================
# 0. Required Packages
# ============================================================================
if (!require("MASS")) install.packages("MASS")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("reshape2")) install.packages("reshape2")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("dplyr")) install.packages("dplyr")

library(MASS)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)

# ============================================================================
# 1. Core Functions: Classical WABA and I-WABA
# ============================================================================

#' Classical WABA Decomposition
#' Decomposes total correlation into Between and Within components
#' using group means E(X|Y)
#'
#' @param X numeric matrix (n x p)
#' @param Y factor or integer vector of group labels (length n)
#' @return list with Total, Between, Within correlation matrices and eta values
classical_waba <- function(X, Y) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  Y <- as.factor(Y)
  groups <- levels(Y)
  K <- length(groups)

  # Grand mean
  grand_mean <- colMeans(X)

  # Total Sum of Squares and Cross Products (centered)
  X_centered <- sweep(X, 2, grand_mean)
  TSS <- crossprod(X_centered)  # t(X_centered) %*% X_centered

  # Between Sum of Squares and Cross Products
  BSS <- matrix(0, p, p)
  for (k in 1:K) {
    idx <- which(Y == groups[k])
    n_k <- length(idx)
    group_mean <- colMeans(X[idx, , drop = FALSE])
    diff <- group_mean - grand_mean
    BSS <- BSS + n_k * outer(diff, diff)
  }

  # Within Sum of Squares and Cross Products
  WSS <- TSS - BSS

  # Convert to correlations
  # Total correlation
  D_total <- sqrt(diag(TSS))
  D_total[D_total == 0] <- 1  # avoid division by zero
  R_total <- TSS / outer(D_total, D_total)

  # Between correlation
  D_between <- sqrt(diag(BSS))
  D_between[D_between == 0] <- 1
  R_between <- BSS / outer(D_between, D_between)

  # Within correlation
  D_within <- sqrt(diag(WSS))
  D_within[D_within == 0] <- 1
  R_within <- WSS / outer(D_within, D_within)

  # Eta values (correlation ratios)
  eta_between <- sqrt(diag(BSS)) / sqrt(diag(TSS))
  eta_within <- sqrt(diag(WSS)) / sqrt(diag(TSS))

  # WABA decomposition: R_total = eta_b * eta_b' * R_between + eta_w * eta_w' * R_within
  R_between_weighted <- outer(eta_between, eta_between) * R_between
  R_within_weighted <- outer(eta_within, eta_within) * R_within

  return(list(
    R_total = R_total,
    R_between = R_between,
    R_within = R_within,
    R_between_weighted = R_between_weighted,
    R_within_weighted = R_within_weighted,
    eta_between = eta_between,
    eta_within = eta_within,
    BSS = BSS,
    WSS = WSS,
    TSS = TSS
  ))
}

#' SDA Covariance (Billard 2008)
#' Computes the covariance between two interval-valued variables
#' using the formula: Cov_SDA = Cov(C, C') + (1/3) Cov(R, R')
#'
#' @param centers_j numeric vector of centers for variable j
#' @param centers_k numeric vector of centers for variable k
#' @param radii_j numeric vector of radii for variable j
#' @param radii_k numeric vector of radii for variable k
#' @return scalar SDA covariance
sda_cov <- function(centers_j, centers_k, radii_j, radii_k) {
  cov_centers <- cov(centers_j, centers_k) * (length(centers_j) - 1) / length(centers_j)
  cov_radii <- cov(radii_j, radii_k) * (length(radii_j) - 1) / length(radii_j)
  return(cov_centers + (1/3) * cov_radii)
}

#' SDA Variance
#' @param centers numeric vector of centers
#' @param radii numeric vector of radii
#' @return scalar SDA variance
sda_var <- function(centers, radii) {
  var_c <- var(centers) * (length(centers) - 1) / length(centers)
  var_r <- var(radii) * (length(radii) - 1) / length(radii)
  return(var_c + (1/3) * var_r)
}

#' Interval-based WABA (I-WABA) Decomposition
#' Uses Int(X|Y) instead of E(X|Y) for the between component
#' Group intervals are defined by (center = group mean, radius = group SD or range/2)
#'
#' @param X numeric matrix (n x p)
#' @param Y factor or integer vector of group labels (length n)
#' @param radius_type character, "sd" for standard deviation, "range" for (max-min)/2
#' @return list with Total, Between (SDA), Within correlation matrices and eta values
iwaba <- function(X, Y, radius_type = "sd") {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  Y <- as.factor(Y)
  groups <- levels(Y)
  K <- length(groups)

  # --- Step 1: Compute group centers and radii ---
  group_centers <- matrix(0, K, p)
  group_radii <- matrix(0, K, p)
  group_sizes <- numeric(K)

  for (k in 1:K) {
    idx <- which(Y == groups[k])
    group_sizes[k] <- length(idx)
    group_centers[k, ] <- colMeans(X[idx, , drop = FALSE])
    if (radius_type == "sd") {
      group_radii[k, ] <- apply(X[idx, , drop = FALSE], 2, sd)
    } else if (radius_type == "range") {
      group_radii[k, ] <- apply(X[idx, , drop = FALSE], 2, function(x) diff(range(x)) / 2)
    }
  }

  # --- Step 2: Compute SDA covariance matrix of Int(X|Y) (Between-SDA) ---
  # This uses weighted versions based on group sizes
  # Weighted centers and radii
  weights <- group_sizes / sum(group_sizes)

  # Weighted mean of centers and radii
  mean_centers <- colSums(weights * group_centers)
  mean_radii <- colSums(weights * group_radii)

  # Between-SDA covariance matrix
  Cov_SDA_between <- matrix(0, p, p)
  for (j in 1:p) {
    for (l in j:p) {
      # Weighted covariance of centers
      cov_c <- sum(weights * (group_centers[, j] - mean_centers[j]) *
                      (group_centers[, l] - mean_centers[l]))
      # Weighted covariance of radii
      cov_r <- sum(weights * (group_radii[, j] - mean_radii[j]) *
                      (group_radii[, l] - mean_radii[l]))
      Cov_SDA_between[j, l] <- cov_c + (1/3) * cov_r
      Cov_SDA_between[l, j] <- Cov_SDA_between[j, l]
    }
  }

  # --- Step 3: Compute Total SDA covariance ---
  # For each observation, treat it as a point interval [x, x] (center=x, radius=0)
  # Then group: Int(X|Y) has center = group mean, radius = group sd
  # Total SDA uses individual observations as intervals
  # Actually, total SDA covariance needs to consider all observations
  # Using the decomposition: Cov_SDA_Total = Cov_SDA_Within + Cov_SDA_Between

  # Within-SDA covariance: average within-group SDA covariance
  # For individual observations within a group, each is a point [x_i, x_i]
  # The SDA covariance of points = standard covariance (radius=0)
  # So Within-SDA = standard Within covariance

  # Standard within-group covariance
  Cov_within <- matrix(0, p, p)
  for (k in 1:K) {
    idx <- which(Y == groups[k])
    X_k <- X[idx, , drop = FALSE]
    X_k_centered <- sweep(X_k, 2, group_centers[k, ])
    Cov_within <- Cov_within + crossprod(X_k_centered)
  }
  Cov_within <- Cov_within / n

  # Total covariance (standard)
  grand_mean <- colMeans(X)
  X_centered <- sweep(X, 2, grand_mean)
  Cov_total_standard <- crossprod(X_centered) / n

  # Between standard covariance (from group means only)
  Cov_between_standard <- Cov_total_standard - Cov_within

  # The total SDA covariance in I-WABA framework
  # Cov_SDA_Total = Cov_within + Cov_SDA_between
  Cov_SDA_total <- Cov_within + Cov_SDA_between

  # --- Step 4: Convert to correlations ---
  # Total SDA correlation
  D_total_sda <- sqrt(diag(Cov_SDA_total))
  D_total_sda[D_total_sda == 0] <- 1
  R_total_sda <- Cov_SDA_total / outer(D_total_sda, D_total_sda)

  # Between SDA correlation
  D_between_sda <- sqrt(diag(Cov_SDA_between))
  D_between_sda[D_between_sda == 0] <- 1
  R_between_sda <- Cov_SDA_between / outer(D_between_sda, D_between_sda)

  # Within correlation (standard)
  D_within <- sqrt(diag(Cov_within))
  D_within[D_within == 0] <- 1
  R_within <- Cov_within / outer(D_within, D_within)

  # Eta values (SDA version)
  eta_between_sda <- sqrt(pmax(diag(Cov_SDA_between), 0)) / sqrt(pmax(diag(Cov_SDA_total), 0))
  eta_within_sda <- sqrt(pmax(diag(Cov_within), 0)) / sqrt(pmax(diag(Cov_SDA_total), 0))

  # Weighted correlation matrices
  R_between_weighted <- outer(eta_between_sda, eta_between_sda) * R_between_sda
  R_within_weighted <- outer(eta_within_sda, eta_within_sda) * R_within

  return(list(
    R_total_sda = R_total_sda,
    R_between_sda = R_between_sda,
    R_within = R_within,
    R_between_weighted = R_between_weighted,
    R_within_weighted = R_within_weighted,
    eta_between_sda = eta_between_sda,
    eta_within_sda = eta_within_sda,
    Cov_SDA_between = Cov_SDA_between,
    Cov_within = Cov_within,
    Cov_SDA_total = Cov_SDA_total,
    Cov_between_standard = Cov_between_standard,
    group_centers = group_centers,
    group_radii = group_radii,
    group_sizes = group_sizes
  ))
}


# ============================================================================
# 2. Data Generation Functions for Simulation Scenarios
# ============================================================================

#' Generate data under heteroscedastic groups with correlated spreads
#' This is the key scenario where I-WABA should outperform classical WABA
#'
#' @param n total sample size
#' @param p number of variables
#' @param K number of groups
#' @param mean_shift magnitude of mean differences across groups
#' @param spread_corr correlation of spreads across variables within groups
#' @param seed random seed
#' @return list with X matrix and Y labels
generate_heteroscedastic <- function(n, p, K, mean_shift = 1.0, spread_corr = 0.8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_per_group <- rep(floor(n / K), K)
  n_per_group[K] <- n - sum(n_per_group[-K])

  # Group means: linearly spaced
  group_means <- matrix(0, K, p)
  for (j in 1:p) {
    group_means[, j] <- seq(-mean_shift, mean_shift, length.out = K) * (1 + 0.2 * (j - 1))
  }

  # Group standard deviations: correlated across variables
  # The key: spreads vary across groups AND are correlated across variables
  base_sd <- seq(0.5, 2.0, length.out = K)  # increasing spread across groups
  group_sds <- matrix(0, K, p)
  for (k in 1:K) {
    # Base spread for this group plus correlated noise
    noise <- mvrnorm(1, mu = rep(0, p),
                     Sigma = spread_corr * matrix(1, p, p) + (1 - spread_corr) * diag(p))
    group_sds[k, ] <- pmax(base_sd[k] + 0.3 * noise, 0.1)
  }

  # Generate data
  X <- matrix(0, n, p)
  Y <- integer(n)
  idx_start <- 1
  for (k in 1:K) {
    idx_end <- idx_start + n_per_group[k] - 1
    for (j in 1:p) {
      X[idx_start:idx_end, j] <- rnorm(n_per_group[k],
                                         mean = group_means[k, j],
                                         sd = group_sds[k, j])
    }
    Y[idx_start:idx_end] <- k
    idx_start <- idx_end + 1
  }

  return(list(X = X, Y = as.factor(Y),
              true_means = group_means, true_sds = group_sds))
}

#' Generate data with spread-only signal (no mean differences)
#' Groups have identical means but different spreads
#' I-WABA should detect the signal, classical WABA should not
#'
#' @param n total sample size
#' @param p number of variables
#' @param K number of groups
#' @param spread_range range of standard deviations across groups
#' @param spread_corr correlation of spreads across variables
#' @param seed random seed
generate_spread_only <- function(n, p, K, spread_range = c(0.5, 3.0),
                                  spread_corr = 0.9, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_per_group <- rep(floor(n / K), K)
  n_per_group[K] <- n - sum(n_per_group[-K])

  # All group means identical (zero)
  group_means <- matrix(0, K, p)

  # Spreads systematically vary and are correlated
  base_sd <- seq(spread_range[1], spread_range[2], length.out = K)
  group_sds <- matrix(0, K, p)
  for (k in 1:K) {
    noise <- mvrnorm(1, mu = rep(0, p),
                     Sigma = spread_corr * matrix(1, p, p) + (1 - spread_corr) * diag(p))
    group_sds[k, ] <- pmax(base_sd[k] + 0.2 * noise, 0.1)
  }

  X <- matrix(0, n, p)
  Y <- integer(n)
  idx_start <- 1
  for (k in 1:K) {
    idx_end <- idx_start + n_per_group[k] - 1
    for (j in 1:p) {
      X[idx_start:idx_end, j] <- rnorm(n_per_group[k],
                                         mean = 0,
                                         sd = group_sds[k, j])
    }
    Y[idx_start:idx_end] <- k
    idx_start <- idx_end + 1
  }

  return(list(X = X, Y = as.factor(Y),
              true_means = group_means, true_sds = group_sds))
}

#' Generate homoscedastic data (control scenario)
#' I-WABA should perform similarly to classical WABA
generate_homoscedastic <- function(n, p, K, mean_shift = 1.0, common_sd = 1.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_per_group <- rep(floor(n / K), K)
  n_per_group[K] <- n - sum(n_per_group[-K])

  group_means <- matrix(0, K, p)
  for (j in 1:p) {
    group_means[, j] <- seq(-mean_shift, mean_shift, length.out = K)
  }

  X <- matrix(0, n, p)
  Y <- integer(n)
  idx_start <- 1
  for (k in 1:K) {
    idx_end <- idx_start + n_per_group[k] - 1
    X[idx_start:idx_end, ] <- mvrnorm(n_per_group[k],
                                       mu = group_means[k, ],
                                       Sigma = common_sd^2 * diag(p))
    Y[idx_start:idx_end] <- k
    idx_start <- idx_end + 1
  }

  return(list(X = X, Y = as.factor(Y),
              true_means = group_means, true_sds = matrix(common_sd, K, p)))
}

#' Generate data with nonlinear variance structure
#' Spread increases with mean (multiplicative heteroscedasticity)
generate_nonlinear_variance <- function(n, p, K, mean_shift = 2.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_per_group <- rep(floor(n / K), K)
  n_per_group[K] <- n - sum(n_per_group[-K])

  group_means <- matrix(0, K, p)
  group_sds <- matrix(0, K, p)
  for (j in 1:p) {
    group_means[, j] <- seq(1, 1 + mean_shift, length.out = K)
    # Spread proportional to mean (coefficient of variation = 0.5)
    group_sds[, j] <- 0.5 * group_means[, j]
  }

  X <- matrix(0, n, p)
  Y <- integer(n)
  idx_start <- 1
  for (k in 1:K) {
    idx_end <- idx_start + n_per_group[k] - 1
    for (j in 1:p) {
      X[idx_start:idx_end, j] <- rnorm(n_per_group[k],
                                         mean = group_means[k, j],
                                         sd = group_sds[k, j])
    }
    Y[idx_start:idx_end] <- k
    idx_start <- idx_end + 1
  }

  return(list(X = X, Y = as.factor(Y),
              true_means = group_means, true_sds = group_sds))
}


# ============================================================================
# 3. Evaluation Metrics
# ============================================================================

#' Compute comparison metrics between Classical WABA and I-WABA
#'
#' @param X data matrix
#' @param Y group labels
#' @return data.frame with comparison metrics
compare_methods <- function(X, Y) {
  waba_result <- classical_waba(X, Y)
  iwaba_result <- iwaba(X, Y, radius_type = "sd")

  p <- ncol(X)

  # Extract upper triangular elements (pairwise correlations)
  upper_idx <- upper.tri(diag(p))

  # Between correlations
  r_between_waba <- waba_result$R_between[upper_idx]
  r_between_iwaba <- iwaba_result$R_between_sda[upper_idx]

  # Between weighted correlations
  r_between_w_waba <- waba_result$R_between_weighted[upper_idx]
  r_between_w_iwaba <- iwaba_result$R_between_weighted[upper_idx]

  # Eta values
  eta_b_waba <- mean(waba_result$eta_between)
  eta_b_iwaba <- mean(iwaba_result$eta_between_sda)

  # Information gain: ratio of I-WABA between variance to WABA between variance
  info_gain <- mean(diag(iwaba_result$Cov_SDA_between)) /
    max(mean(diag(iwaba_result$Cov_between_standard)), 1e-10)

  # Radius contribution: proportion of between-SDA variance due to radius
  radius_var <- mean(pmax(diag(iwaba_result$Cov_SDA_between) -
                            diag(iwaba_result$Cov_between_standard), 0))
  total_between_sda_var <- mean(diag(iwaba_result$Cov_SDA_between))
  radius_proportion <- ifelse(total_between_sda_var > 0,
                               radius_var / total_between_sda_var, 0)

  data.frame(
    mean_abs_between_waba = mean(abs(r_between_waba)),
    mean_abs_between_iwaba = mean(abs(r_between_iwaba)),
    eta_between_waba = eta_b_waba,
    eta_between_iwaba = eta_b_iwaba,
    information_gain = info_gain,
    radius_proportion = radius_proportion
  )
}


# ============================================================================
# 4. Main Simulation Study
# ============================================================================

run_simulation <- function(n_reps = 100,
                            K_values = c(2, 3, 5, 10),
                            n_values = c(100, 200, 500, 1000),
                            p_values = c(2, 5, 10, 50),
                            scenarios = c("heteroscedastic", "spread_only",
                                          "homoscedastic", "nonlinear_variance"),
                            verbose = TRUE) {

  results <- list()
  counter <- 0

  total_combos <- length(scenarios) * length(K_values) * length(n_values) * length(p_values)

  for (scenario in scenarios) {
    for (K in K_values) {
      for (n in n_values) {
        # Ensure n >= K * 5 (at least 5 per group)
        if (n < K * 5) next

        for (p in p_values) {
          counter <- counter + 1
          if (verbose) cat(sprintf("[%d/%d] Scenario=%s, K=%d, n=%d, p=%d\n",
                                   counter, total_combos, scenario, K, n, p))

          rep_results <- vector("list", n_reps)

          for (rep in 1:n_reps) {
            # Generate data
            dat <- switch(scenario,
                          "heteroscedastic" = generate_heteroscedastic(n, p, K, seed = rep * 1000 + counter),
                          "spread_only" = generate_spread_only(n, p, K, seed = rep * 1000 + counter),
                          "homoscedastic" = generate_homoscedastic(n, p, K, seed = rep * 1000 + counter),
                          "nonlinear_variance" = generate_nonlinear_variance(n, p, K, seed = rep * 1000 + counter))

            # Compare methods
            tryCatch({
              metrics <- compare_methods(dat$X, dat$Y)
              metrics$scenario <- scenario
              metrics$K <- K
              metrics$n <- n
              metrics$p <- p
              metrics$rep <- rep
              rep_results[[rep]] <- metrics
            }, error = function(e) {
              if (verbose) cat(sprintf("  Error in rep %d: %s\n", rep, e$message))
            })
          }

          results <- c(results, rep_results[!sapply(rep_results, is.null)])
        }
      }
    }
  }

  results_df <- do.call(rbind, results)
  return(results_df)
}


# ============================================================================
# 5. Run Simulations
# ============================================================================

cat("============================================================\n")
cat("Starting I-WABA Simulation Study\n")
cat("============================================================\n\n")

# Run with moderate settings for tractability
sim_results <- run_simulation(
  n_reps = 100,
  K_values = c(2, 3, 5, 10),
  n_values = c(100, 200, 500, 1000),
  p_values = c(2, 5, 10, 50),
  scenarios = c("heteroscedastic", "spread_only", "homoscedastic", "nonlinear_variance"),
  verbose = TRUE
)

# Save results
saveRDS(sim_results, file = "simulation_results.rds")
cat("\nSimulation results saved to simulation_results.rds\n")


# ============================================================================
# 6. Summary Statistics
# ============================================================================

cat("\n============================================================\n")
cat("Summary Statistics by Scenario and K\n")
cat("============================================================\n\n")

summary_stats <- sim_results %>%
  group_by(scenario, K) %>%
  summarize(
    mean_info_gain = mean(information_gain, na.rm = TRUE),
    sd_info_gain = sd(information_gain, na.rm = TRUE),
    mean_radius_prop = mean(radius_proportion, na.rm = TRUE),
    mean_eta_waba = mean(eta_between_waba, na.rm = TRUE),
    mean_eta_iwaba = mean(eta_between_iwaba, na.rm = TRUE),
    eta_improvement = mean(eta_between_iwaba - eta_between_waba, na.rm = TRUE),
    .groups = "drop"
  )

print(as.data.frame(summary_stats))

# Summary by scenario, K, and n
summary_by_n <- sim_results %>%
  group_by(scenario, K, n) %>%
  summarize(
    mean_info_gain = mean(information_gain, na.rm = TRUE),
    mean_radius_prop = mean(radius_proportion, na.rm = TRUE),
    eta_improvement = mean(eta_between_iwaba - eta_between_waba, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================================
# 7. Visualization
# ============================================================================

# --- Figure 1: Information Gain vs K by Scenario ---
fig1 <- ggplot(summary_stats, aes(x = factor(K), y = mean_info_gain,
                                    fill = scenario)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_info_gain - sd_info_gain,
                     ymax = mean_info_gain + sd_info_gain),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.5) +
  labs(title = "Information Gain of I-WABA over Classical WABA",
       subtitle = "Values > 1 indicate I-WABA captures more between-group information",
       x = "Number of Classes (K)",
       y = "Information Gain (I-WABA / WABA)",
       fill = "Scenario") +
  theme_bw(base_size = 12) +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "bottom")

ggsave("fig_info_gain_vs_K.pdf", fig1, width = 10, height = 6)

# --- Figure 2: Radius Contribution vs K ---
fig2 <- ggplot(summary_stats, aes(x = factor(K), y = mean_radius_prop,
                                    color = scenario, group = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(title = "Radius Contribution to Between-Group Variance",
       subtitle = "Proportion of between-SDA variance attributable to group spread differences",
       x = "Number of Classes (K)",
       y = "Radius Proportion",
       color = "Scenario") +
  theme_bw(base_size = 12) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "bottom")

ggsave("fig_radius_contribution.pdf", fig2, width = 10, height = 6)

# --- Figure 3: Eta Improvement by K and n ---
fig3 <- ggplot(summary_by_n %>% filter(scenario %in% c("heteroscedastic", "spread_only")),
               aes(x = factor(K), y = eta_improvement, color = factor(n), group = factor(n))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~scenario, scales = "free_y") +
  labs(title = "Improvement in Eta (Correlation Ratio) by I-WABA",
       subtitle = expression(eta[between]^{I-WABA} - eta[between]^{WABA}),
       x = "Number of Classes (K)",
       y = "Eta Improvement",
       color = "Sample Size (n)") +
  theme_bw(base_size = 12) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom")

ggsave("fig_eta_improvement.pdf", fig3, width = 12, height = 6)

# --- Figure 4: Heatmap of Information Gain (K vs p) ---
heatmap_data <- sim_results %>%
  filter(scenario == "heteroscedastic") %>%
  group_by(K, p) %>%
  summarize(mean_info_gain = mean(information_gain, na.rm = TRUE), .groups = "drop")

fig4 <- ggplot(heatmap_data, aes(x = factor(K), y = factor(p), fill = mean_info_gain)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", mean_info_gain)), size = 4) +
  scale_fill_gradient2(low = "white", mid = "steelblue", high = "darkblue",
                        midpoint = median(heatmap_data$mean_info_gain)) +
  labs(title = "Information Gain: I-WABA vs WABA (Heteroscedastic Scenario)",
       x = "Number of Classes (K)",
       y = "Number of Variables (p)",
       fill = "Info Gain") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

ggsave("fig_heatmap_info_gain.pdf", fig4, width = 8, height = 6)

# --- Figure 5: Boxplot comparison across scenarios ---
fig5_data <- sim_results %>%
  filter(n == 500, p == 10) %>%
  select(scenario, K, information_gain) %>%
  mutate(K = factor(K))

fig5 <- ggplot(fig5_data, aes(x = K, y = information_gain, fill = scenario)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = "Distribution of Information Gain (n=500, p=10)",
       x = "Number of Classes (K)",
       y = "Information Gain",
       fill = "Scenario") +
  theme_bw(base_size = 12) +
  scale_fill_brewer(palette = "Pastel1") +
  theme(legend.position = "bottom")

ggsave("fig_boxplot_info_gain.pdf", fig5, width = 10, height = 6)

cat("\n============================================================\n")
cat("All simulation figures saved.\n")
cat("============================================================\n")


# ============================================================================
# 8. Generate LaTeX Table
# ============================================================================

cat("\n============================================================\n")
cat("LaTeX Table: Summary of Simulation Results\n")
cat("============================================================\n\n")

table_data <- sim_results %>%
  filter(n == 500, p == 10) %>%
  group_by(scenario, K) %>%
  summarize(
    info_gain_mean = sprintf("%.3f", mean(information_gain, na.rm = TRUE)),
    info_gain_sd = sprintf("%.3f", sd(information_gain, na.rm = TRUE)),
    radius_prop = sprintf("%.3f", mean(radius_proportion, na.rm = TRUE)),
    eta_waba = sprintf("%.3f", mean(eta_between_waba, na.rm = TRUE)),
    eta_iwaba = sprintf("%.3f", mean(eta_between_iwaba, na.rm = TRUE)),
    .groups = "drop"
  )

cat("\\begin{table}[ht]\n")
cat("\\centering\n")
cat("\\caption{Simulation results: I-WABA vs Classical WABA ($n=500$, $p=10$, 100 replications)}\n")
cat("\\label{tab:sim_results}\n")
cat("\\begin{tabular}{llccccc}\n")
cat("\\toprule\n")
cat("Scenario & $K$ & Info Gain (Mean) & Info Gain (SD) & Radius Prop. & $\\bar{\\eta}_B^{\\text{WABA}}$ & $\\bar{\\eta}_B^{\\text{I-WABA}}$ \\\\\n")
cat("\\midrule\n")

for (i in 1:nrow(table_data)) {
  row <- table_data[i, ]
  cat(sprintf("%s & %s & %s & %s & %s & %s & %s \\\\\n",
              row$scenario, row$K, row$info_gain_mean, row$info_gain_sd,
              row$radius_prop, row$eta_waba, row$eta_iwaba))
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")

cat("\n============================================================\n")
cat("Simulation study complete.\n")
cat("============================================================\n")
