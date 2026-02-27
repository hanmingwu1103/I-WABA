################################################################################
## I-WABA: Core Functions
## Title: Interval-based Within-And-Between Analysis
## Author: Han-Ming Wu (National Chengchi University)
## Reference: Wu, H.-M. (2026). I-WABA: Interval-based Within-And-Between
##            Analysis. Statistical Analysis and Data Mining.
################################################################################

#' Classical WABA Eta Correlation Ratios
#'
#' Computes between-group (eta_B) and within-group (eta_W) correlation ratios
#' for a single variable using classical WABA (group means).
#'
#' @param x numeric vector of observations
#' @param g factor or vector of group labels
#' @return list with eta_B, eta_W, E (ratio), BSS, WSS, TSS
classical_waba_eta <- function(x, g) {
  grand_mean <- mean(x)
  groups <- split(x, g)
  n_total <- length(x)

  BSS <- sum(sapply(groups, function(grp) {
    length(grp) * (mean(grp) - grand_mean)^2
  }))
  WSS <- sum(sapply(groups, function(grp) {
    sum((grp - mean(grp))^2)
  }))
  TSS <- BSS + WSS

  eta_B <- sqrt(BSS / TSS)
  eta_W <- sqrt(WSS / TSS)
  E <- eta_B / eta_W

  list(eta_B = eta_B, eta_W = eta_W, E = E,
       BSS = BSS, WSS = WSS, TSS = TSS)
}

#' SDA Between-Group Covariance Matrix
#'
#' Computes the Symbolic Data Analysis (SDA) between-group covariance:
#'   Cov_SDA = Cov(Centers) + (1/3) * Cov(Radii)
#' where Centers = weighted group means, Radii = weighted group SDs.
#'
#' @param X numeric matrix (n x p) of observations
#' @param g factor or vector of group labels
#' @return list with cov_sda, cov_centers, cov_radii, centers, radii, weights
sda_between_covariance <- function(X, g) {
  X <- as.matrix(X)
  groups <- split(data.frame(X), g)
  K <- length(groups)
  p <- ncol(X)
  n <- nrow(X)

  centers <- matrix(0, K, p)
  radii <- matrix(0, K, p)
  weights <- numeric(K)
  group_names <- names(groups)

  for (k in seq_len(K)) {
    grp <- as.matrix(groups[[k]])
    n_k <- nrow(grp)
    centers[k, ] <- colMeans(grp)
    radii[k, ] <- apply(grp, 2, sd) * sqrt((n_k - 1) / n_k)  # population SD
    weights[k] <- n_k
  }

  weights_norm <- weights / sum(weights)

  # Weighted covariance of centers
  center_mean <- colSums(centers * weights_norm)
  center_dev <- sweep(centers, 2, center_mean)
  cov_centers <- t(center_dev) %*% diag(weights_norm) %*% center_dev

  # Weighted covariance of radii
  radii_mean <- colSums(radii * weights_norm)
  radii_dev <- sweep(radii, 2, radii_mean)
  cov_radii <- t(radii_dev) %*% diag(weights_norm) %*% radii_dev

  # SDA covariance
  cov_sda <- cov_centers + (1/3) * cov_radii

  rownames(centers) <- group_names
  rownames(radii) <- group_names

  list(cov_sda = cov_sda, cov_centers = cov_centers, cov_radii = cov_radii,
       centers = centers, radii = radii, weights = weights,
       group_names = group_names)
}

#' Classical Between-Group Covariance Matrix
#'
#' @param X numeric matrix (n x p)
#' @param g group labels
#' @return between-group covariance matrix (p x p)
classical_between_covariance <- function(X, g) {
  X <- as.matrix(X)
  grand_mean <- colMeans(X)
  groups <- split(data.frame(X), g)
  n <- nrow(X)
  p <- ncol(X)

  cov_between <- matrix(0, p, p)
  for (grp in groups) {
    grp <- as.matrix(grp)
    n_k <- nrow(grp)
    dev <- colMeans(grp) - grand_mean
    cov_between <- cov_between + (n_k / n) * outer(dev, dev)
  }
  cov_between
}

#' I-WABA Eta Correlation Ratios
#'
#' Computes I-WABA eta_B and eta_W for a single variable using
#' SDA between-group covariance.
#'
#' @param x numeric vector
#' @param g group labels
#' @param sda_between_var scalar: SDA between-group variance for this variable
#' @return list with eta_B_sda, eta_W_sda, E_sda
iwaba_eta <- function(x, g, sda_between_var) {
  waba <- classical_waba_eta(x, g)
  n <- length(x)

  BSS_SDA <- n * sda_between_var
  WSS <- waba$WSS
  TSS_SDA <- BSS_SDA + WSS

  eta_B_sda <- sqrt(BSS_SDA / TSS_SDA)
  eta_W_sda <- sqrt(WSS / TSS_SDA)
  E_sda <- eta_B_sda / eta_W_sda

  list(eta_B_sda = eta_B_sda, eta_W_sda = eta_W_sda, E_sda = E_sda)
}

#' Full I-WABA Analysis
#'
#' Performs the complete I-WABA I (Entity Analysis) and computes
#' information gain, radius proportion, and per-variable eta comparison.
#'
#' @param X numeric matrix (n x p) or data.frame
#' @param g group labels (factor or vector)
#' @return list with detailed results
iwaba_analysis <- function(X, g) {
  X <- as.matrix(X)
  p <- ncol(X)
  n <- nrow(X)
  K <- length(unique(g))

  # SDA covariance
  sda <- sda_between_covariance(X, g)
  cov_classical <- classical_between_covariance(X, g)

  # Information gain
  info_gain <- sum(diag(sda$cov_sda)) / sum(diag(cov_classical))
  radius_prop <- (1/3) * sum(diag(sda$cov_radii)) / sum(diag(sda$cov_sda))

  # Per-variable eta analysis
  results <- data.frame(
    variable = if (!is.null(colnames(X))) colnames(X) else paste0("V", 1:p),
    eta_B_waba = numeric(p),
    eta_W_waba = numeric(p),
    E_waba = numeric(p),
    eta_B_iwaba = numeric(p),
    eta_W_iwaba = numeric(p),
    E_iwaba = numeric(p),
    delta_eta = numeric(p),
    E_ratio = numeric(p),
    stringsAsFactors = FALSE
  )

  for (j in 1:p) {
    w <- classical_waba_eta(X[, j], g)
    iw <- iwaba_eta(X[, j], g, sda$cov_sda[j, j])

    results$eta_B_waba[j] <- w$eta_B
    results$eta_W_waba[j] <- w$eta_W
    results$E_waba[j] <- w$E
    results$eta_B_iwaba[j] <- iw$eta_B_sda
    results$eta_W_iwaba[j] <- iw$eta_W_sda
    results$E_iwaba[j] <- iw$E_sda
    results$delta_eta[j] <- iw$eta_B_sda - w$eta_B
    results$E_ratio[j] <- iw$E_sda / w$E
  }

  list(
    n = n, K = K, p = p,
    info_gain = info_gain,
    radius_prop = radius_prop,
    eta_results = results,
    sda = sda,
    cov_classical = cov_classical
  )
}

#' WABA I Classification
#'
#' Classifies variables as Between, Equivocal, or Within based on E-test ratio.
#' Thresholds follow Dansereau et al. (1984): 15-degree rule.
#'
#' @param E numeric E-test ratio
#' @return character classification
classify_e <- function(E) {
  ifelse(E >= 1.30, "Between",
         ifelse(E <= 1/1.30, "Within", "Equivocal"))
}

#' Print I-WABA Analysis Summary
#'
#' @param result output from iwaba_analysis()
print_iwaba_summary <- function(result) {
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("I-WABA Analysis Summary\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat(sprintf("  n = %d, K = %d, p = %d\n", result$n, result$K, result$p))
  cat(sprintf("  Information Gain G = %.4f\n", result$info_gain))
  cat(sprintf("  Radius Proportion  = %.4f (%.1f%%)\n",
              result$radius_prop, result$radius_prop * 100))
  cat(sprintf("  Mean η_B (WABA):   %.4f\n", mean(result$eta_results$eta_B_waba)))
  cat(sprintf("  Mean η_B (I-WABA): %.4f\n", mean(result$eta_results$eta_B_iwaba)))
  cat(sprintf("  Mean Δη_B:         %.4f\n", mean(result$eta_results$delta_eta)))
  cat(sprintf("  Max Δη_B:          %.4f (%s)\n",
              max(result$eta_results$delta_eta),
              result$eta_results$variable[which.max(result$eta_results$delta_eta)]))
  cat("───────────────────────────────────────────────────────────────\n")
}
