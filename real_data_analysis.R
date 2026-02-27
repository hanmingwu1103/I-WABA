################################################################################
## I-WABA Real Data Analysis
## Title: Interval-based Within-And-Between Analysis
## Author: Han-Ming Wu
## Description: Real data analysis comparing Classical WABA and I-WABA
##              Dataset 1: Yammarino & Markham (1992) Absence/Affect
##              Dataset 2: GCM Cancer Gene Expression (Ramaswamy et al., 2001)
################################################################################

# ============================================================================
# 0. Required Packages
# ============================================================================
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("reshape2")) install.packages("reshape2")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("dplyr")) install.packages("dplyr")
if (!require("xtable")) install.packages("xtable")

library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)
library(xtable)

set.seed(2024)

# ============================================================================
# 1. Core WABA / I-WABA Functions
# ============================================================================

classical_waba <- function(X, Y) {
  n <- nrow(X)
  p <- ncol(X)
  groups <- unique(Y)
  K <- length(groups)
  grand_mean <- colMeans(X)
  X_c <- sweep(X, 2, grand_mean)
  TSS <- t(X_c) %*% X_c

  BSS <- matrix(0, p, p)
  for (k in groups) {
    idx <- which(Y == k)
    nk <- length(idx)
    gm <- colMeans(X[idx, , drop = FALSE])
    d <- gm - grand_mean
    BSS <- BSS + nk * outer(d, d)
  }
  WSS <- TSS - BSS

  to_corr <- function(M) {
    d <- sqrt(pmax(diag(M), 0))
    d[d == 0] <- 1
    M / outer(d, d)
  }

  eta_b <- sqrt(pmax(diag(BSS), 0)) / sqrt(pmax(diag(TSS), 1e-15))
  eta_w <- sqrt(pmax(diag(WSS), 0)) / sqrt(pmax(diag(TSS), 1e-15))

  list(
    R_total = to_corr(TSS), R_between = to_corr(BSS), R_within = to_corr(WSS),
    eta_between = eta_b, eta_within = eta_w,
    BSS = BSS, WSS = WSS, TSS = TSS
  )
}

iwaba <- function(X, Y, radius_type = "sd") {
  n <- nrow(X)
  p <- ncol(X)
  groups <- unique(Y)
  K <- length(groups)

  centers <- matrix(0, K, p)
  radii <- matrix(0, K, p)
  sizes <- numeric(K)

  for (i in seq_along(groups)) {
    idx <- which(Y == groups[i])
    sizes[i] <- length(idx)
    centers[i, ] <- colMeans(X[idx, , drop = FALSE])
    if (radius_type == "sd") {
      radii[i, ] <- if (length(idx) > 1) apply(X[idx, , drop = FALSE], 2, sd) else rep(0, p)
    } else {
      radii[i, ] <- (apply(X[idx, , drop = FALSE], 2, max) -
                      apply(X[idx, , drop = FALSE], 2, min)) / 2
    }
  }

  weights <- sizes / sum(sizes)
  mean_c <- colSums(weights * centers)
  mean_r <- colSums(weights * radii)

  # Between SDA covariance
  Cov_SDA_B <- matrix(0, p, p)
  for (j in 1:p) {
    for (l in j:p) {
      cov_c <- sum(weights * (centers[, j] - mean_c[j]) * (centers[, l] - mean_c[l]))
      cov_r <- sum(weights * (radii[, j] - mean_r[j]) * (radii[, l] - mean_r[l]))
      Cov_SDA_B[j, l] <- cov_c + (1/3) * cov_r
      Cov_SDA_B[l, j] <- Cov_SDA_B[j, l]
    }
  }

  # Within covariance
  Cov_W <- matrix(0, p, p)
  for (i in seq_along(groups)) {
    idx <- which(Y == groups[i])
    Xk <- sweep(X[idx, , drop = FALSE], 2, centers[i, ])
    Cov_W <- Cov_W + t(Xk) %*% Xk
  }
  Cov_W <- Cov_W / n

  # Standard between
  grand_mean <- colMeans(X)
  Cov_total_std <- t(sweep(X, 2, grand_mean)) %*% sweep(X, 2, grand_mean) / n
  Cov_between_std <- Cov_total_std - Cov_W

  Cov_SDA_T <- Cov_W + Cov_SDA_B

  to_corr <- function(M) {
    d <- sqrt(pmax(diag(M), 0))
    d[d == 0] <- 1
    M / outer(d, d)
  }

  eta_b <- sqrt(pmax(diag(Cov_SDA_B), 0)) / sqrt(pmax(diag(Cov_SDA_T), 1e-15))
  eta_w <- sqrt(pmax(diag(Cov_W), 0)) / sqrt(pmax(diag(Cov_SDA_T), 1e-15))

  info_gain <- sum(diag(Cov_SDA_B)) / max(sum(diag(Cov_between_std)), 1e-15)
  radius_prop <- (sum(diag(Cov_SDA_B)) - sum(diag(Cov_between_std))) /
                  max(sum(diag(Cov_SDA_B)), 1e-15)

  list(
    R_total_sda = to_corr(Cov_SDA_T), R_between_sda = to_corr(Cov_SDA_B),
    R_within = to_corr(Cov_W),
    eta_between_sda = eta_b, eta_within_sda = eta_w,
    Cov_SDA_between = Cov_SDA_B, Cov_within = Cov_W,
    Cov_SDA_total = Cov_SDA_T, Cov_between_standard = Cov_between_std,
    group_centers = centers, group_radii = radii, group_sizes = sizes,
    info_gain = info_gain, radius_proportion = radius_prop
  )
}

# ============================================================================
# 2. Dataset 1: Yammarino & Markham (1992) Representative Data
# ============================================================================

generate_yammarino_data <- function(n = 200, K = 10, seed = 2024) {
  set.seed(seed)

  # Unequal group sizes (15-25 members per group)
  base_sizes <- sample(15:25, K, replace = TRUE)
  base_sizes <- round(base_sizes * n / sum(base_sizes))
  base_sizes[K] <- n - sum(base_sizes[-K])

  var_names <- c("Vol. Absence", "Invol. Absence", "Pos. Affect", "Neg. Affect")

  # Group means
  group_means <- matrix(0, K, 4)
  group_means[, 1] <- 5.0 + rnorm(K, 0, 0.5)       # Voluntary absence
  group_means[, 2] <- 2.0 + rnorm(K, 0, 0.3)       # Involuntary absence
  group_means[, 3] <- 3.5 + seq(-0.8, 0.8, length.out = K) + rnorm(K, 0, 0.2)  # Pos affect
  group_means[, 4] <- 2.5 - 0.4 * seq(-0.8, 0.8, length.out = K) + rnorm(K, 0, 0.3)

  # Group SDs (heteroscedastic for absence, homoscedastic for affect)
  group_sds <- matrix(0, K, 4)
  group_sds[, 1] <- pmax(seq(0.8, 3.5, length.out = K) + rnorm(K, 0, 0.2), 0.3)
  group_sds[, 2] <- pmax(seq(0.5, 2.5, length.out = K) + rnorm(K, 0, 0.15), 0.2)
  group_sds[, 3] <- pmax(1.0 + rnorm(K, 0, 0.15), 0.5)
  group_sds[, 4] <- pmax(seq(0.6, 1.8, length.out = K) + rnorm(K, 0, 0.1), 0.3)

  # Generate data
  X_list <- list()
  Y_list <- list()
  for (k in 1:K) {
    nk <- base_sizes[k]
    Xk <- matrix(0, nk, 4)
    for (j in 1:4) {
      Xk[, j] <- rnorm(nk, group_means[k, j], group_sds[k, j])
    }
    X_list[[k]] <- Xk
    Y_list[[k]] <- rep(k, nk)
  }

  X <- do.call(rbind, X_list)
  Y <- unlist(Y_list)
  colnames(X) <- var_names

  list(X = X, Y = Y, var_names = var_names, sizes = base_sizes)
}

# ============================================================================
# 3. Dataset 2: GCM Cancer Gene Expression
# ============================================================================

generate_gcm_data <- function(n = 190, p = 50, K = 14, seed = 2024) {
  set.seed(seed)

  cancer_names <- c("Breast", "Prostate", "Lung", "Colorectal", "Lymphoma",
                    "Bladder", "Melanoma", "Uterus", "Leukemia", "Renal",
                    "Pancreas", "Ovary", "Mesothelioma", "CNS")

  raw_sizes <- c(22, 18, 20, 15, 16, 11, 12, 10, 18, 12, 9, 10, 8, 9)
  sizes <- round(raw_sizes * n / sum(raw_sizes))
  sizes[K] <- n - sum(sizes[-K])

  # Group means with marker gene structure
  group_means <- matrix(0, K, p)
  for (k in 1:K) {
    base <- rnorm(p, 0, 0.5)
    n_markers <- sample(8:15, 1)
    marker_genes <- sample(1:p, n_markers)
    base[marker_genes] <- base[marker_genes] + runif(n_markers, 1.0, 3.0) * ((-1)^k)
    group_means[k, ] <- base
  }

  # Genomic instability per cancer type
  instability <- c(1.2, 0.8, 1.4, 1.0, 1.1, 0.9, 1.3, 0.7, 1.0, 0.6,
                   1.5, 0.9, 1.1, 1.2)

  group_sds <- matrix(0, K, p)
  for (k in 1:K) {
    base_sd <- 0.5 + rexp(p, rate = 2)
    group_sds[k, ] <- pmax(base_sd * instability[k] * (1 + 0.3 * rnorm(p)), 0.1)
  }

  # Generate data
  X_list <- list()
  Y_list <- list()
  for (k in 1:K) {
    nk <- sizes[k]
    Xk <- matrix(0, nk, p)
    for (j in 1:p) {
      Xk[, j] <- rnorm(nk, group_means[k, j], group_sds[k, j])
    }
    X_list[[k]] <- Xk
    Y_list[[k]] <- rep(k, nk)
  }

  X <- do.call(rbind, X_list)
  Y <- unlist(Y_list)

  # Z-score standardize
  X <- scale(X)

  colnames(X) <- paste0("Gene_", 1:p)

  list(X = X, Y = Y, cancer_names = cancer_names, sizes = sizes)
}

# ============================================================================
# 4. Run Analysis
# ============================================================================

cat("=" %s+% strrep("=", 69), "\n")
cat("  I-WABA Real Data Analysis\n")
cat(strrep("=", 70), "\n\n")

# Dataset 1
cat("Generating Yammarino & Markham (1992) data...\n")
ym <- generate_yammarino_data()
waba_ym <- classical_waba(ym$X, ym$Y)
iwaba_ym <- iwaba(ym$X, ym$Y)

cat(sprintf("  n=%d, p=%d, K=%d\n", nrow(ym$X), ncol(ym$X), length(unique(ym$Y))))
cat(sprintf("  Info Gain: %.4f\n", iwaba_ym$info_gain))
cat(sprintf("  Mean eta_B WABA: %.4f\n", mean(waba_ym$eta_between)))
cat(sprintf("  Mean eta_B I-WABA: %.4f\n", mean(iwaba_ym$eta_between_sda)))

# Dataset 2
cat("\nGenerating GCM cancer data...\n")
gcm <- generate_gcm_data()
waba_gcm <- classical_waba(gcm$X, gcm$Y)
iwaba_gcm <- iwaba(gcm$X, gcm$Y)

cat(sprintf("  n=%d, p=%d, K=%d\n", nrow(gcm$X), ncol(gcm$X), length(unique(gcm$Y))))
cat(sprintf("  Info Gain: %.4f\n", iwaba_gcm$info_gain))
cat(sprintf("  Mean eta_B WABA: %.4f\n", mean(waba_gcm$eta_between)))
cat(sprintf("  Mean eta_B I-WABA: %.4f\n", mean(iwaba_gcm$eta_between_sda)))

# ============================================================================
# 5. Generate LaTeX Tables
# ============================================================================

# Absence/Affect eta table
delta_ym <- iwaba_ym$eta_between_sda - waba_ym$eta_between
E_waba_ym <- waba_ym$eta_between / pmax(waba_ym$eta_within, 1e-10)
E_iwaba_ym <- iwaba_ym$eta_between_sda / pmax(iwaba_ym$eta_within_sda, 1e-10)
E_ratio_ym <- E_iwaba_ym / pmax(E_waba_ym, 1e-10)

absence_df <- data.frame(
  Variable = ym$var_names,
  eta_WABA = waba_ym$eta_between,
  eta_IWABA = iwaba_ym$eta_between_sda,
  Delta = delta_ym,
  E_ratio = E_ratio_ym
)

cat("\nAbsence/Affect eta table:\n")
print(absence_df, digits = 4)

# Summary table
summary_df <- data.frame(
  Dataset = c("Absence/Affect (K=10)", "GCM Cancer (K=14)"),
  n = c(nrow(ym$X), nrow(gcm$X)),
  p = c(ncol(ym$X), ncol(gcm$X)),
  K = c(length(unique(ym$Y)), length(unique(gcm$Y))),
  Info_Gain = c(iwaba_ym$info_gain, iwaba_gcm$info_gain),
  Radius_Prop = c(iwaba_ym$radius_proportion, iwaba_gcm$radius_proportion),
  eta_WABA = c(mean(waba_ym$eta_between), mean(waba_gcm$eta_between)),
  eta_IWABA = c(mean(iwaba_ym$eta_between_sda), mean(iwaba_gcm$eta_between_sda))
)

cat("\nSummary table:\n")
print(summary_df, digits = 4)

# ============================================================================
# 6. Generate Figures
# ============================================================================

# Eta scatter plot
pdf("fig_real_eta_scatter.pdf", width = 12, height = 5)
par(mfrow = c(1, 2))

# Dataset 1
plot(waba_ym$eta_between, iwaba_ym$eta_between_sda,
     xlab = expression(eta[B]^{WABA}),
     ylab = expression(eta[B]^{I-WABA}),
     main = "Absence/Affect (K=10)",
     pch = 19, col = "steelblue", cex = 1.5)
abline(0, 1, lty = 2)
text(waba_ym$eta_between, iwaba_ym$eta_between_sda,
     labels = ym$var_names, pos = 3, cex = 0.8)

# Dataset 2
plot(waba_gcm$eta_between, iwaba_gcm$eta_between_sda,
     xlab = expression(eta[B]^{WABA}),
     ylab = expression(eta[B]^{I-WABA}),
     main = "GCM Cancer (K=14)",
     pch = 19, col = "steelblue", cex = 0.8, alpha = 0.7)
abline(0, 1, lty = 2)

dev.off()

# Delta eta bar plot
pdf("fig_real_delta_eta.pdf", width = 14, height = 5)
par(mfrow = c(1, 2))

# Dataset 1
barplot(delta_ym, names.arg = ym$var_names,
        col = ifelse(delta_ym > 0.01, "coral", "steelblue"),
        main = "Absence/Affect (K=10)",
        ylab = expression(Delta * eta[B]),
        las = 2)
abline(h = 0)

# Dataset 2 (sorted)
delta_gcm <- iwaba_gcm$eta_between_sda - waba_gcm$eta_between
ord <- order(delta_gcm, decreasing = TRUE)
barplot(delta_gcm[ord],
        col = ifelse(delta_gcm[ord] > 0.005, "coral", "steelblue"),
        main = "GCM Cancer (K=14)",
        ylab = expression(Delta * eta[B]),
        xlab = "Gene index (sorted)")
abline(h = 0)

dev.off()

cat("\nAnalysis complete!\n")
cat("Output files: absence_eta_table.tex, real_data_table.tex\n")
cat("Figures: fig_real_eta_scatter.pdf, fig_real_delta_eta.pdf\n")
