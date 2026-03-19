################################################################################
## I-WABA: Core Functions (v2 — Full Level I–V Analysis)
## Title: Interval-based Within-And-Between Analysis
## Author: Han-Ming Wu (National Chengchi University)
## Reference: Wu, H.-M. (2026). I-WABA: Interval-based Within-And-Between
##            Analysis. Statistical Analysis and Data Mining.
##
## Key change from v1:
##   Uses Billard's endpoint covariance (Eq. 21):
##     V_R^(j) = (1/3) * sum_k w_k * s_kj^2   [uncentered second moment]
##   instead of the centered-radius version:
##     (1/3) * Cov(R_j, R_j')
##
## Implements all five I-WABA levels:
##   Level I   — Entity Analysis (eta, E-test, info gain)
##   Level II  — Functional Relationship (weighted correlation decomposition)
##   Level III — Consistency Check (Level I vs Level II agreement)
##   Level IV  — Dispersion-Source Diagnostics (center vs radius attribution)
##   Level V   — Dispersion Association (radius correlations across variables)
################################################################################

# =============================================================================
# Utility: safe division
# =============================================================================
safe_div <- function(a, b, fallback = 0) {
  ifelse(abs(b) > 1e-15, a / b, fallback)
}

# =============================================================================
# Utility: covariance-to-correlation
# =============================================================================
cov_to_corr <- function(M) {
  d <- sqrt(pmax(diag(M), 0))
  d[d == 0] <- 1e-15
  M / outer(d, d)
}

# =============================================================================
# E-test classification (15-degree rule, Dansereau et al. 1984)
# =============================================================================
classify_e <- function(E) {
  ifelse(E >= 1.30, "Between",
         ifelse(E <= 1 / 1.30, "Within", "Equivocal"))
}


# =============================================================================
# Level I: Entity Analysis
# =============================================================================
#' Compute group-level quantities: centers, radii, weights
#'
#' @param X numeric matrix (n x p)
#' @param g factor or vector of group labels (length n)
#' @return list with centers (K x p), radii (K x p), n_k, w_k, K, p, n
group_quantities <- function(X, g) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  g <- as.factor(g)
  groups <- levels(g)
  K <- length(groups)

  centers <- matrix(0, K, p)
  radii   <- matrix(0, K, p)    # population SD: s_kj
  n_k     <- numeric(K)

  for (ki in seq_len(K)) {
    idx <- which(g == groups[ki])
    n_k[ki] <- length(idx)
    X_k <- X[idx, , drop = FALSE]
    centers[ki, ] <- colMeans(X_k)
    # Population SD (ddof = 0)
    if (n_k[ki] > 1) {
      radii[ki, ] <- sqrt(colMeans(sweep(X_k, 2, centers[ki, ])^2))
    }
  }

  w_k <- n_k / n
  rownames(centers) <- groups
  rownames(radii)   <- groups

  list(centers = centers, radii = radii, n_k = n_k, w_k = w_k,
       K = K, p = p, n = n, groups = groups)
}


#' I-WABA Level I Entity Analysis
#'
#' Computes V_C, V_R (Billard), V_W, eta ratios, E-test, info gain.
#'
#' @param X numeric matrix (n x p) or data.frame
#' @param g group labels
#' @return list with Level I results
iwaba_level1 <- function(X, g) {
  X <- as.matrix(X)
  gq <- group_quantities(X, g)
  n <- gq$n; p <- gq$p; K <- gq$K
  w_k <- gq$w_k; centers <- gq$centers; radii <- gq$radii

  # Grand weighted mean of centers
  center_mean <- colSums(w_k * centers)

  # V_C^(j) = sum_k w_k (C_kj - C_bar_j)^2
  V_C <- colSums(w_k * sweep(centers, 2, center_mean)^2)

  # V_R^(j) = (1/3) * sum_k w_k * s_kj^2  [Billard's formulation, Eq. 21]
  V_R <- (1/3) * colSums(w_k * radii^2)

  # V_W^(j) = pooled within-group variance (W_jj / n)
  g_factor <- as.factor(g)
  V_W <- numeric(p)
  for (ki in seq_len(K)) {
    idx <- which(g_factor == gq$groups[ki])
    X_k <- X[idx, , drop = FALSE]
    V_W <- V_W + colSums(sweep(X_k, 2, centers[ki, ])^2) / n
  }

  # I-WABA between and augmented total
  Var_between_int <- V_C + V_R
  Var_aug <- V_C + V_R + V_W

  # I-WABA eta ratios
  eta_B_int  <- sqrt(safe_div(Var_between_int, Var_aug))
  eta_W_int  <- sqrt(safe_div(V_W, Var_aug))

  # Classical WABA eta
  Var_total_classical <- V_C + V_W
  eta_B_waba <- sqrt(safe_div(V_C, Var_total_classical))
  eta_W_waba <- sqrt(safe_div(V_W, Var_total_classical))

  # E-test ratios
  E_waba <- sqrt(safe_div(V_C, V_W))
  E_int  <- sqrt(safe_div(V_C + V_R, V_W))

  # Information gain G (trace-based)
  info_gain   <- safe_div(sum(Var_between_int), sum(V_C), fallback = Inf)
  radius_prop <- safe_div(sum(V_R), sum(Var_between_int))

  # Delta eta
  delta_eta <- eta_B_int - eta_B_waba

  # E-test classifications
  class_waba  <- classify_e(E_waba)
  class_iwaba <- classify_e(E_int)
  n_reclassified <- sum(class_waba != class_iwaba)

  # Per-variable results data.frame
  varnames <- if (!is.null(colnames(X))) colnames(X) else paste0("V", seq_len(p))
  eta_results <- data.frame(
    variable    = varnames,
    eta_B_waba  = eta_B_waba,
    eta_W_waba  = eta_W_waba,
    E_waba      = E_waba,
    eta_B_iwaba = eta_B_int,
    eta_W_iwaba = eta_W_int,
    E_iwaba     = E_int,
    delta_eta   = delta_eta,
    E_ratio     = safe_div(E_int, E_waba),
    class_waba  = class_waba,
    class_iwaba = class_iwaba,
    stringsAsFactors = FALSE
  )

  list(
    n = n, K = K, p = p,
    V_C = V_C, V_R = V_R, V_W = V_W,
    info_gain = info_gain,
    radius_prop = radius_prop,
    eta_results = eta_results,
    gq = gq
  )
}


# =============================================================================
# Level IV: Dispersion-Source Diagnostics
# =============================================================================
#' @param gq group_quantities output
#' @param V_C, V_R, V_W from Level I
#' @return list with pi_R, V_R_het, E_R_het, etc.
iwaba_level4 <- function(gq, V_C, V_R, V_W) {
  K <- gq$K; p <- gq$p; w_k <- gq$w_k; radii <- gq$radii

  Var_between_int <- V_C + V_R

  # pi_R: radius share per variable
  pi_R <- safe_div(V_R, Var_between_int)

  # E_C^(Int) and E_R^(Int)
  E_C_int <- sqrt(safe_div(V_C, V_W))
  E_R_int <- sqrt(safe_div(V_R, V_W))

  # V_{R,het}^(j) = (1/3) * sum_k w_k (s_kj - s_bar_j)^2
  s_bar <- colSums(w_k * radii)   # weighted mean radius
  V_R_het <- (1/3) * colSums(w_k * sweep(radii, 2, s_bar)^2)

  E_R_het_int <- sqrt(safe_div(V_R_het, V_W))

  list(
    pi_R = pi_R,
    V_R_het = V_R_het,
    E_C_int = E_C_int,
    E_R_int = E_R_int,
    E_R_het_int = E_R_het_int,
    s_bar = s_bar
  )
}


# =============================================================================
# Level II: Functional Relationship Analysis
# =============================================================================
#' @param X data matrix (n x p)
#' @param g group labels
#' @param gq group_quantities
#' @param eta_B_int, eta_W_int I-WABA eta vectors (length p)
#' @param eta_B_waba, eta_W_waba classical WABA eta vectors (length p)
#' @return list with covariance matrices and weighted components
iwaba_level2 <- function(X, g, gq, eta_B_int, eta_W_int, eta_B_waba, eta_W_waba) {
  X <- as.matrix(X)
  n <- gq$n; p <- gq$p; K <- gq$K
  w_k <- gq$w_k; centers <- gq$centers; radii <- gq$radii
  g_factor <- as.factor(g)

  center_mean <- colSums(w_k * centers)
  C_tilde <- sweep(centers, 2, center_mean)   # centered group means

  # Between-group interval covariance (Billard, Eq. 21)
  Cov_between_int <- matrix(0, p, p)
  radius_product  <- matrix(0, p, p)
  for (ki in seq_len(K)) {
    Cov_between_int <- Cov_between_int + w_k[ki] * tcrossprod(C_tilde[ki, ])
    radius_product  <- radius_product  + w_k[ki] * tcrossprod(radii[ki, ])
  }
  Cov_between_int <- Cov_between_int + radius_product / 3

  # Classical between-group covariance (centers only)
  Cov_between_classical <- matrix(0, p, p)
  for (ki in seq_len(K)) {
    Cov_between_classical <- Cov_between_classical + w_k[ki] * tcrossprod(C_tilde[ki, ])
  }

  # Within-group covariance
  Cov_within <- matrix(0, p, p)
  for (ki in seq_len(K)) {
    idx <- which(g_factor == gq$groups[ki])
    dev <- sweep(X[idx, , drop = FALSE], 2, centers[ki, ])
    Cov_within <- Cov_within + crossprod(dev) / n
  }

  # Correlations
  r_between_int      <- cov_to_corr(Cov_between_int)
  r_between_classical <- cov_to_corr(Cov_between_classical)
  r_within           <- cov_to_corr(Cov_within)

  # Pairwise weighted components
  if (p >= 2) {
    pairs <- combn(seq_len(p), 2)
    n_pairs <- ncol(pairs)

    bw_int <- numeric(n_pairs)
    ww_int <- numeric(n_pairs)
    bw_waba <- numeric(n_pairs)
    ww_waba <- numeric(n_pairs)

    for (pi_idx in seq_len(n_pairs)) {
      i <- pairs[1, pi_idx]
      j <- pairs[2, pi_idx]

      bw_int[pi_idx]  <- eta_B_int[i]  * eta_B_int[j]  * r_between_int[i, j]
      ww_int[pi_idx]  <- eta_W_int[i]  * eta_W_int[j]  * r_within[i, j]
      bw_waba[pi_idx] <- eta_B_waba[i] * eta_B_waba[j] * r_between_classical[i, j]
      ww_waba[pi_idx] <- eta_W_waba[i] * eta_W_waba[j] * r_within[i, j]
    }

    mean_abs_bw_int  <- mean(abs(bw_int))
    mean_abs_ww_int  <- mean(abs(ww_int))
    mean_abs_bw_waba <- mean(abs(bw_waba))
    mean_abs_ww_waba <- mean(abs(ww_waba))
  } else {
    pairs <- matrix(nrow = 2, ncol = 0)
    n_pairs <- 0
    bw_int <- ww_int <- bw_waba <- ww_waba <- numeric(0)
    mean_abs_bw_int <- mean_abs_ww_int <- 0
    mean_abs_bw_waba <- mean_abs_ww_waba <- 0
  }

  list(
    Cov_between_int = Cov_between_int,
    Cov_between_classical = Cov_between_classical,
    Cov_within = Cov_within,
    r_between_int = r_between_int,
    r_between_classical = r_between_classical,
    r_within = r_within,
    pairs = pairs, n_pairs = n_pairs,
    bw_int = bw_int, ww_int = ww_int,
    bw_waba = bw_waba, ww_waba = ww_waba,
    mean_abs_bw_int = mean_abs_bw_int,
    mean_abs_ww_int = mean_abs_ww_int,
    mean_abs_bw_waba = mean_abs_bw_waba,
    mean_abs_ww_waba = mean_abs_ww_waba
  )
}


# =============================================================================
# Level III: Consistency Check
# =============================================================================
#' @param class_iwaba, class_waba character vectors of E-test classifications
#' @param bw_int, ww_int, bw_waba, ww_waba pairwise weighted components
#' @param pairs 2 x n_pairs matrix of pair indices
#' @return list with consistency_rate_int and consistency_rate_waba
iwaba_level3 <- function(class_iwaba, class_waba,
                          bw_int, ww_int, bw_waba, ww_waba, pairs) {
  n_pairs <- ncol(pairs)
  if (n_pairs == 0) return(list(consistency_rate_int = 1, consistency_rate_waba = 1))

  consistency_int  <- 0
  consistency_waba <- 0

  for (pi_idx in seq_len(n_pairs)) {
    i <- pairs[1, pi_idx]
    j <- pairs[2, pi_idx]

    both_B_int <- (class_iwaba[i] == "Between") & (class_iwaba[j] == "Between")
    l2_B_int   <- abs(bw_int[pi_idx]) > abs(ww_int[pi_idx])
    if (both_B_int == l2_B_int) consistency_int <- consistency_int + 1

    both_B_waba <- (class_waba[i] == "Between") & (class_waba[j] == "Between")
    l2_B_waba   <- abs(bw_waba[pi_idx]) > abs(ww_waba[pi_idx])
    if (both_B_waba == l2_B_waba) consistency_waba <- consistency_waba + 1
  }

  list(
    consistency_rate_int  = consistency_int / n_pairs,
    consistency_rate_waba = consistency_waba / n_pairs
  )
}


# =============================================================================
# Level V: Dispersion Association
# =============================================================================
#' @param gq group_quantities
#' @param s_bar weighted mean radii (from Level IV)
#' @return list with r_R_billard and r_R_het per pair, plus diagnostic metadata
iwaba_level5 <- function(gq, s_bar) {
  K <- gq$K; p <- gq$p; w_k <- gq$w_k; radii <- gq$radii

  s_tilde <- sweep(radii, 2, s_bar)  # centered radii
  is_k2_degenerate <- (K == 2)

  if (p < 2) {
    return(list(
      r_R_billard = numeric(0),
      r_R_het = numeric(0),
      r_R_billard_defined = logical(0),
      r_R_het_defined = logical(0),
      r_R_het_sign_only = logical(0),
      is_k2_degenerate = is_k2_degenerate,
      n_pairs = 0,
      n_r_R_billard_defined = 0,
      n_r_R_het_defined = 0,
      n_r_R_het_undefined = 0,
      prop_r_R_het_positive = NA_real_,
      prop_r_R_het_negative = NA_real_,
      prop_r_R_het_undefined = NA_real_,
      mean_r_R_billard = NA_real_,
      mean_r_R_het_signed = NA_real_,
      mean_abs_r_R_billard = NA_real_,
      mean_abs_r_R_het = NA_real_
    ))
  }

  pairs <- combn(seq_len(p), 2)
  n_pairs <- ncol(pairs)

  r_R_billard <- rep(NA_real_, n_pairs)
  r_R_het     <- rep(NA_real_, n_pairs)
  r_R_billard_defined <- rep(FALSE, n_pairs)
  r_R_het_defined     <- rep(FALSE, n_pairs)

  for (pi_idx in seq_len(n_pairs)) {
    i <- pairs[1, pi_idx]
    j <- pairs[2, pi_idx]

    # Billard radius association (magnitude-inclusive)
    num_b   <- sum(w_k * radii[, i] * radii[, j])
    denom_b <- sqrt(sum(w_k * radii[, i]^2) * sum(w_k * radii[, j]^2))
    if (denom_b > 1e-15) {
      r_R_billard[pi_idx] <- min(max(num_b / denom_b, 0), 1)
      r_R_billard_defined[pi_idx] <- TRUE
    }

    # Heterogeneity-based radius correlation (centered)
    num_h   <- sum(w_k * s_tilde[, i] * s_tilde[, j])
    denom_h <- sqrt(sum(w_k * s_tilde[, i]^2) * sum(w_k * s_tilde[, j]^2))
    if (denom_h > 1e-15) {
      r_R_het[pi_idx] <- max(min(num_h / denom_h, 1), -1)
      r_R_het_defined[pi_idx] <- TRUE
    }
  }

  n_r_R_billard_defined <- sum(r_R_billard_defined)
  n_r_R_het_defined <- sum(r_R_het_defined)
  n_r_R_het_undefined <- n_pairs - n_r_R_het_defined

  r_R_het_defined_values <- r_R_het[r_R_het_defined]
  if (n_r_R_het_defined > 0) {
    prop_r_R_het_positive <- mean(r_R_het_defined_values > 0)
    prop_r_R_het_negative <- mean(r_R_het_defined_values < 0)
  } else {
    prop_r_R_het_positive <- NA_real_
    prop_r_R_het_negative <- NA_real_
  }
  prop_r_R_het_undefined <- n_r_R_het_undefined / n_pairs

  if (n_r_R_billard_defined > 0) {
    mean_r_R_billard <- mean(r_R_billard[r_R_billard_defined])
  } else {
    mean_r_R_billard <- NA_real_
  }
  if (!is_k2_degenerate && n_r_R_het_defined > 0) {
    mean_r_R_het_signed <- mean(r_R_het_defined_values)
  } else {
    mean_r_R_het_signed <- NA_real_
  }

  list(
    r_R_billard = r_R_billard,
    r_R_het = r_R_het,
    r_R_billard_defined = r_R_billard_defined,
    r_R_het_defined = r_R_het_defined,
    r_R_het_sign_only = rep(is_k2_degenerate, n_pairs) & r_R_het_defined,
    is_k2_degenerate = is_k2_degenerate,
    n_pairs = n_pairs,
    n_r_R_billard_defined = n_r_R_billard_defined,
    n_r_R_het_defined = n_r_R_het_defined,
    n_r_R_het_undefined = n_r_R_het_undefined,
    prop_r_R_het_positive = prop_r_R_het_positive,
    prop_r_R_het_negative = prop_r_R_het_negative,
    prop_r_R_het_undefined = prop_r_R_het_undefined,
    mean_r_R_billard = mean_r_R_billard,
    mean_r_R_het_signed = mean_r_R_het_signed,
    mean_abs_r_R_billard = mean(abs(r_R_billard), na.rm = TRUE),
    mean_abs_r_R_het = mean(abs(r_R_het), na.rm = TRUE)
  )
}


# =============================================================================
# Full I-WABA Analysis (all 5 levels)
# =============================================================================
#' Performs the complete I-WABA Level I–V analysis.
#'
#' @param X numeric matrix (n x p) or data.frame
#' @param g group labels (factor or vector)
#' @return list with results from all five levels
iwaba_full <- function(X, g) {
  X <- as.matrix(X)

  # Level I
  l1 <- iwaba_level1(X, g)

  # Level IV
  l4 <- iwaba_level4(l1$gq, l1$V_C, l1$V_R, l1$V_W)

  # Level II
  eta_B_int  <- l1$eta_results$eta_B_iwaba
  eta_W_int  <- l1$eta_results$eta_W_iwaba
  eta_B_waba <- l1$eta_results$eta_B_waba
  eta_W_waba <- l1$eta_results$eta_W_waba

  l2 <- iwaba_level2(X, g, l1$gq, eta_B_int, eta_W_int, eta_B_waba, eta_W_waba)

  # Level III
  l3 <- iwaba_level3(
    l1$eta_results$class_iwaba, l1$eta_results$class_waba,
    l2$bw_int, l2$ww_int, l2$bw_waba, l2$ww_waba, l2$pairs
  )

  # Level V
  l5 <- iwaba_level5(l1$gq, l4$s_bar)

  list(
    level1 = l1,
    level2 = l2,
    level3 = l3,
    level4 = l4,
    level5 = l5,
    gq = l1$gq
  )
}


# =============================================================================
# Legacy wrapper: iwaba_analysis (for backward compatibility)
# =============================================================================
iwaba_analysis <- function(X, g) {
  result <- iwaba_full(X, g)
  l1 <- result$level1

  list(
    n = l1$n, K = l1$K, p = l1$p,
    info_gain = l1$info_gain,
    radius_prop = l1$radius_prop,
    eta_results = l1$eta_results,
    sda = list(
      centers = l1$gq$centers,
      radii   = l1$gq$radii,
      weights = l1$gq$n_k,
      group_names = l1$gq$groups
    ),
    level2 = result$level2,
    level3 = result$level3,
    level4 = result$level4,
    level5 = result$level5
  )
}


# =============================================================================
# Classical WABA (for comparison)
# =============================================================================
classical_waba <- function(X, g) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  g <- as.factor(g)
  groups <- levels(g)
  K <- length(groups)
  grand_mean <- colMeans(X)

  X_centered <- sweep(X, 2, grand_mean)
  TSS <- crossprod(X_centered)

  BSS <- matrix(0, p, p)
  for (ki in seq_len(K)) {
    idx <- which(g == groups[ki])
    n_k <- length(idx)
    gm  <- colMeans(X[idx, , drop = FALSE])
    d   <- gm - grand_mean
    BSS <- BSS + n_k * tcrossprod(d)
  }
  WSS <- TSS - BSS

  eta_between <- sqrt(pmax(diag(BSS), 0)) / sqrt(pmax(diag(TSS), 1e-15))
  eta_within  <- sqrt(pmax(diag(WSS), 0)) / sqrt(pmax(diag(TSS), 1e-15))

  R_between <- cov_to_corr(BSS)
  R_within  <- cov_to_corr(WSS)

  R_between_weighted <- outer(eta_between, eta_between) * R_between
  R_within_weighted  <- outer(eta_within, eta_within)   * R_within

  list(
    R_total = cov_to_corr(TSS),
    R_between = R_between,
    R_within = R_within,
    R_between_weighted = R_between_weighted,
    R_within_weighted = R_within_weighted,
    eta_between = eta_between,
    eta_within = eta_within,
    BSS = BSS, WSS = WSS, TSS = TSS
  )
}


# =============================================================================
# Print Summary
# =============================================================================
print_iwaba_summary <- function(result) {
  # Accept either iwaba_analysis output or iwaba_full output
  if (!is.null(result$level1)) {
    l1 <- result$level1
    eta <- l1$eta_results
    info_gain <- l1$info_gain
    radius_prop <- l1$radius_prop
    n <- l1$n; K <- l1$K; p <- l1$p
  } else {
    eta <- result$eta_results
    info_gain <- result$info_gain
    radius_prop <- result$radius_prop
    n <- result$n; K <- result$K; p <- result$p
  }

  cat("================================================================\n")
  cat("I-WABA Analysis Summary (Billard's endpoint covariance)\n")
  cat("================================================================\n")
  cat(sprintf("  n = %d, K = %d, p = %d\n", n, K, p))
  cat(sprintf("  Information Gain G = %.4f\n", info_gain))
  cat(sprintf("  Radius Proportion  = %.4f (%.1f%%)\n",
              radius_prop, radius_prop * 100))
  cat(sprintf("  Mean eta_B (WABA):   %.4f\n", mean(eta$eta_B_waba)))
  cat(sprintf("  Mean eta_B (I-WABA): %.4f\n", mean(eta$eta_B_iwaba)))
  cat(sprintf("  Mean Delta_eta_B:    %.4f\n", mean(eta$delta_eta)))
  cat(sprintf("  Max  Delta_eta_B:    %.4f (%s)\n",
              max(eta$delta_eta),
              eta$variable[which.max(eta$delta_eta)]))
  cat(sprintf("  Reclassified:        %d / %d\n",
              sum(eta$class_waba != eta$class_iwaba), p))
  cat("----------------------------------------------------------------\n")
}
