################################################################################
## I-WABA Applied Inference Helpers
##
## Shared helpers for real-data analyses:
##   - dominance contrasts
##   - within-group bootstrap wrappers
##   - Brown-Forsythe / Levene-type tests
##   - BH/FDR adjustment
##   - global heterogeneity screening summary
##
## This file does not modify the core I-WABA formulas. It assumes
## `iwaba_functions.R` has already been sourced.
################################################################################

iwaba_pair_labels <- function(pairs, varnames) {
  if (is.null(dim(pairs)) || ncol(pairs) == 0) {
    return(character(0))
  }
  apply(pairs, 2, function(idx) {
    paste0(varnames[idx[1]], "--", varnames[idx[2]])
  })
}

iwaba_entity_dominance_contrast <- function(result_or_level1, j = NULL) {
  l1 <- if (!is.null(result_or_level1$level1)) result_or_level1$level1 else result_or_level1
  delta <- l1$V_C + l1$V_R - l1$V_W

  if (!is.null(l1$eta_results$variable)) {
    names(delta) <- l1$eta_results$variable
  }

  if (is.null(j)) {
    delta
  } else {
    delta[j]
  }
}

iwaba_pairwise_dominance_contrast <- function(result_or_level2, pairs = NULL, varnames = NULL) {
  l2 <- if (!is.null(result_or_level2$level2)) result_or_level2$level2 else result_or_level2
  delta <- l2$bw_int^2 - l2$ww_int^2

  if (is.null(pairs)) {
    pairs <- l2$pairs
  }
  if (!is.null(varnames) && length(delta) > 0) {
    names(delta) <- iwaba_pair_labels(pairs, varnames)
  }

  delta
}

iwaba_ci_excludes_value <- function(lower, upper, value = 0) {
  if (is.na(lower) || is.na(upper)) {
    NA
  } else {
    value < lower || value > upper
  }
}

iwaba_bh_adjust <- function(p_values, alpha = 0.05) {
  p_adj <- p.adjust(p_values, method = "BH")
  data.frame(
    p_value = p_values,
    p_adj = p_adj,
    reject = !is.na(p_adj) & p_adj <= alpha,
    alpha = alpha,
    stringsAsFactors = FALSE
  )
}

iwaba_global_t_het <- function(result_or_level4) {
  l4 <- if (!is.null(result_or_level4$level4)) result_or_level4$level4 else result_or_level4
  sum(l4$V_R_het, na.rm = TRUE)
}

iwaba_brown_forsythe_test <- function(x, g, center = c("median", "mean")) {
  center <- match.arg(center)
  g <- as.factor(g)

  if (length(x) != length(g)) {
    stop("x and g must have the same length.")
  }

  loc_fn <- if (center == "median") stats::median else mean
  group_loc <- tapply(x, g, loc_fn)
  deviations <- abs(x - group_loc[g])
  group_sizes <- as.numeric(table(g))
  group_means <- as.numeric(tapply(deviations, g, mean))
  grand_mean <- mean(deviations)

  ss_between <- sum(group_sizes * (group_means - grand_mean)^2)
  ss_within <- sum(tapply(deviations, g, function(z) sum((z - mean(z))^2)))

  df1 <- nlevels(g) - 1L
  df2 <- length(x) - nlevels(g)

  if (df1 <= 0L || df2 <= 0L || ss_within <= 0) {
    return(list(
      statistic = NA_real_,
      df1 = df1,
      df2 = df2,
      p_value = NA_real_,
      center = center
    ))
  }

  statistic <- (ss_between / df1) / (ss_within / df2)
  p_value <- 1 - pf(statistic, df1, df2)

  list(
    statistic = statistic,
    df1 = df1,
    df2 = df2,
    p_value = p_value,
    center = center
  )
}

iwaba_brown_forsythe_matrix <- function(X, g, varnames = NULL, center = c("median", "mean"), alpha = 0.05) {
  X <- as.matrix(X)
  center <- match.arg(center)
  p <- ncol(X)

  if (is.null(varnames)) {
    varnames <- if (!is.null(colnames(X))) colnames(X) else paste0("V", seq_len(p))
  }

  rows <- vector("list", p)
  for (j in seq_len(p)) {
    bf <- iwaba_brown_forsythe_test(X[, j], g, center = center)
    rows[[j]] <- data.frame(
      variable = varnames[j],
      statistic = bf$statistic,
      df1 = bf$df1,
      df2 = bf$df2,
      p_value = bf$p_value,
      center = bf$center,
      stringsAsFactors = FALSE
    )
  }

  results <- do.call(rbind, rows)
  bh <- iwaba_bh_adjust(results$p_value, alpha = alpha)
  cbind(results, bh[, c("p_adj", "reject"), drop = FALSE])
}

iwaba_bootstrap_interval <- function(samples, conf_level = 0.95, min_finite = 50L) {
  finite <- is.finite(samples)
  n_finite <- sum(finite)

  if (n_finite < min_finite) {
    return(list(
      lower = NA_real_,
      upper = NA_real_,
      n_finite = n_finite,
      defined = FALSE
    ))
  }

  probs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
  ci <- stats::quantile(samples[finite], probs = probs, na.rm = TRUE, names = FALSE, type = 7)

  list(
    lower = unname(ci[1]),
    upper = unname(ci[2]),
    n_finite = n_finite,
    defined = TRUE
  )
}

iwaba_within_group_bootstrap_vector <- function(X, g, statistic_fn,
                                                B = 999L,
                                                conf_level = 0.95,
                                                seed = NULL,
                                                min_finite = 50L,
                                                null_value = NULL,
                                                progress = FALSE) {
  X <- as.matrix(X)
  g <- as.factor(g)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  observed <- statistic_fn(X, g)
  stat_names <- names(observed)
  observed <- as.numeric(observed)
  if (is.null(stat_names)) {
    stat_names <- paste0("stat_", seq_along(observed))
  }
  names(observed) <- stat_names

  group_indices <- split(seq_len(nrow(X)), g)
  boot_matrix <- matrix(NA_real_, nrow = B, ncol = length(observed))
  colnames(boot_matrix) <- stat_names

  pb <- NULL
  if (isTRUE(progress)) {
    pb <- txtProgressBar(min = 0, max = B, style = 3)
  }

  for (b in seq_len(B)) {
    boot_idx <- unlist(lapply(group_indices, function(idx) {
      sample(idx, length(idx), replace = TRUE)
    }), use.names = FALSE)

    boot_stats <- statistic_fn(X[boot_idx, , drop = FALSE], g[boot_idx])
    boot_stats <- as.numeric(boot_stats)
    if (length(boot_stats) != length(observed)) {
      stop("Bootstrap statistic length changed across resamples.")
    }
    boot_matrix[b, ] <- boot_stats

    if (!is.null(pb)) {
      setTxtProgressBar(pb, b)
    }
  }

  if (!is.null(pb)) {
    close(pb)
  }

  ci_rows <- vector("list", length(observed))
  for (j in seq_along(observed)) {
    ci <- iwaba_bootstrap_interval(boot_matrix[, j], conf_level = conf_level, min_finite = min_finite)
    ci_rows[[j]] <- data.frame(
      stat_name = stat_names[j],
      estimate = observed[j],
      ci_lower = ci$lower,
      ci_upper = ci$upper,
      ci_defined = ci$defined,
      n_finite = ci$n_finite,
      boot_defined_rate = ci$n_finite / B,
      reject_null = if (is.null(null_value)) {
        NA
      } else {
        iwaba_ci_excludes_value(ci$lower, ci$upper, value = null_value)
      },
      stringsAsFactors = FALSE
    )
  }

  list(
    observed = observed,
    boot_matrix = boot_matrix,
    ci = do.call(rbind, ci_rows),
    B = B,
    conf_level = conf_level
  )
}

iwaba_within_group_bootstrap_scalar <- function(X, g, statistic_fn,
                                                B = 999L,
                                                conf_level = 0.95,
                                                seed = NULL,
                                                min_finite = 50L,
                                                null_value = NULL,
                                                progress = FALSE) {
  wrapped_fn <- function(Xb, gb) {
    value <- statistic_fn(Xb, gb)
    if (length(value) != 1L) {
      stop("Scalar bootstrap wrapper requires statistic_fn to return length 1.")
    }
    stats::setNames(as.numeric(value), "stat")
  }

  result <- iwaba_within_group_bootstrap_vector(
    X = X,
    g = g,
    statistic_fn = wrapped_fn,
    B = B,
    conf_level = conf_level,
    seed = seed,
    min_finite = min_finite,
    null_value = null_value,
    progress = progress
  )

  list(
    estimate = result$observed[1],
    ci = result$ci[1, , drop = FALSE],
    boot_samples = result$boot_matrix[, 1]
  )
}

iwaba_rhet_bootstrap <- function(X, g, pair, varnames = NULL,
                                 B = 999L,
                                 conf_level = 0.95,
                                 seed = NULL,
                                 min_finite = 50L,
                                 progress = FALSE) {
  X <- as.matrix(X)
  if (length(pair) != 2L) {
    stop("pair must contain exactly two variable indices.")
  }

  if (is.null(varnames)) {
    varnames <- if (!is.null(colnames(X))) colnames(X) else paste0("V", seq_len(ncol(X)))
  }
  pair_label <- paste0(varnames[pair[1]], "--", varnames[pair[2]])

  stat_fn <- function(Xb, gb) {
    result <- iwaba_full(Xb, gb)
    pair_idx <- which(result$level2$pairs[1, ] == pair[1] & result$level2$pairs[2, ] == pair[2])
    if (length(pair_idx) != 1L) {
      stop("Could not locate requested pair in bootstrap output.")
    }
    stats::setNames(result$level5$r_R_het[pair_idx], pair_label)
  }

  boot <- iwaba_within_group_bootstrap_vector(
    X = X,
    g = g,
    statistic_fn = stat_fn,
    B = B,
    conf_level = conf_level,
    seed = seed,
    min_finite = min_finite,
    null_value = 0,
    progress = progress
  )

  cbind(
    data.frame(pair_label = pair_label, stringsAsFactors = FALSE),
    boot$ci
  )
}
