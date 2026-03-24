################################################################################
## I-WABA Real Data Analysis: Bliese-Halverson (bh1996) Dataset
## Staged bh1996 rerun for the real-data milestone
################################################################################

source("iwaba_functions.R")
source("iwaba_inference_helpers.R")

parse_arg_value <- function(args, prefix, default) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) == 0L) {
    default
  } else {
    sub(paste0("^", prefix), "", hit[1])
  }
}

timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
}

append_runlog_entry <- function(runlog_file,
                                subtask,
                                files_changed,
                                commands_run,
                                outputs_created,
                                seeds_used,
                                what_remains) {
  con <- file(runlog_file, open = "a", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  writeLines(c(
    "",
    sprintf("## Milestone: bh1996 / %s", subtask),
    "",
    "- Files changed"
  ), con)
  writeLines(paste0("  - `", files_changed, "`"), con)

  writeLines(c("", "- Commands run"), con)
  writeLines(paste0("  - `", commands_run, "`"), con)

  writeLines(c("", "- Outputs created"), con)
  writeLines(paste0("  - `", outputs_created, "`"), con)

  writeLines(c("", "- Seeds used"), con)
  writeLines(paste0("  - `", seeds_used, "`"), con)

  writeLines(c("", "- What remains"), con)
  writeLines(paste0("  - ", what_remains), con)
}

bind_named_rows <- function(dfs) {
  dfs <- Filter(Negate(is.null), dfs)
  dfs <- Filter(function(x) nrow(x) > 0L, dfs)
  if (length(dfs) == 0L) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  all_cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  aligned <- lapply(dfs, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    if (length(missing_cols) > 0L) {
      for (col_name in missing_cols) {
        df[[col_name]] <- NA
      }
    }
    df[, all_cols, drop = FALSE]
  })

  do.call(rbind, aligned)
}

save_dual_output <- function(csv_object,
                             rds_object,
                             csv_filename,
                             output_dir = "results",
                             row.names = FALSE) {
  csv_path <- iwaba_output_path(csv_filename, output_dir = output_dir, create_dir = TRUE)
  rds_path <- iwaba_output_path(
    sub("\\.csv$", ".rds", csv_filename, ignore.case = TRUE),
    output_dir = output_dir,
    create_dir = TRUE
  )

  write.csv(csv_object, csv_path, row.names = row.names)
  saveRDS(rds_object, rds_path)

  list(csv = csv_path, rds = rds_path, timestamp = timestamp_now())
}

register_dual_output <- function(output_registry, label, paths) {
  rbind(
    output_registry,
    data.frame(
      label = label,
      type = c("csv", "rds"),
      path = c(paths$csv, paths$rds),
      timestamp = paths$timestamp,
      stringsAsFactors = FALSE
    )
  )
}

register_pdf_output <- function(output_registry, label, pdf_path) {
  rbind(
    output_registry,
    data.frame(
      label = label,
      type = "pdf",
      path = pdf_path,
      timestamp = timestamp_now(),
      stringsAsFactors = FALSE
    )
  )
}

format_output_registry <- function(output_registry) {
  if (nrow(output_registry) == 0L) {
    return("")
  }

  paste(
    sprintf(
      "%s [%s] @ %s",
      basename(output_registry$path),
      output_registry$type,
      output_registry$timestamp
    ),
    collapse = "; "
  )
}

build_group_size_summary <- function(g) {
  n_k <- as.numeric(table(g))
  q <- stats::quantile(n_k, probs = c(0.25, 0.50, 0.75), names = FALSE, type = 7)

  data.frame(
    dataset = "bh1996",
    grouping_variable = "GRP",
    K = nlevels(g),
    n = length(g),
    min_n_k = min(n_k),
    q1_n_k = unname(q[1]),
    median_n_k = unname(q[2]),
    mean_n_k = mean(n_k),
    q3_n_k = unname(q[3]),
    max_n_k = max(n_k),
    stringsAsFactors = FALSE
  )
}

build_entity_df <- function(result, X, g, variables) {
  l1 <- result$level1
  eta <- l1$eta_results
  delta_E <- iwaba_entity_dominance_contrast(result)
  icc_values <- vapply(variables, function(v) {
    waba <- classical_waba(X[, v, drop = FALSE], g)
    unname(diag(waba$BSS) / diag(waba$TSS))
  }, numeric(1))

  data.frame(
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
}

build_pair_df <- function(result, variables) {
  l2 <- result$level2
  pair_labels <- iwaba_pair_labels(l2$pairs, variables)
  delta_A <- iwaba_pairwise_dominance_contrast(result, varnames = variables)

  data.frame(
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
}

build_level3_df <- function(result, pair_df) {
  data.frame(
    consistency_rate_waba = result$level3$consistency_rate_waba,
    consistency_rate_iwaba = result$level3$consistency_rate_int,
    n_pairs = nrow(pair_df),
    n_dominance_flips = sum(pair_df$flipped),
    stringsAsFactors = FALSE
  )
}

build_level4_df <- function(result, variables) {
  l4 <- result$level4

  data.frame(
    variable = variables,
    pi_R = l4$pi_R,
    E_C_int = l4$E_C_int,
    E_R_int = l4$E_R_int,
    E_R_het_int = l4$E_R_het_int,
    V_R_het = l4$V_R_het,
    stringsAsFactors = FALSE
  )
}

build_level5_df <- function(result, variables) {
  l2 <- result$level2
  l5 <- result$level5
  pair_labels <- iwaba_pair_labels(l2$pairs, variables)

  data.frame(
    pair_label = pair_labels,
    var1 = variables[l2$pairs[1, ]],
    var2 = variables[l2$pairs[2, ]],
    r_R_billard = l5$r_R_billard,
    r_R_billard_defined = l5$r_R_billard_defined,
    r_R_het = l5$r_R_het,
    r_R_het_defined = l5$r_R_het_defined,
    stringsAsFactors = FALSE
  )
}

build_summary_df <- function(result, entity_df, pair_df, level3_df, level4_df, branch_label) {
  l1 <- result$level1
  l5 <- result$level5

  data.frame(
    branch = branch_label,
    n = l1$n,
    K = l1$K,
    p = l1$p,
    G = l1$info_gain,
    radius_prop = l1$radius_prop,
    mean_eta_B_waba = mean(entity_df$eta_B_waba),
    mean_eta_B_int = mean(entity_df$eta_B_int),
    mean_delta_eta_B = mean(entity_df$delta_eta_B),
    reclassified = sum(entity_df$class_waba != entity_df$class_iwaba),
    mean_pi_R = mean(level4_df$pi_R),
    mean_E_R_het_int = mean(level4_df$E_R_het_int),
    mean_r_R_billard = l5$mean_r_R_billard,
    mean_r_R_het_signed = l5$mean_r_R_het_signed,
    T_het = iwaba_global_t_het(result),
    level3_consistency_waba = level3_df$consistency_rate_waba,
    level3_consistency_iwaba = level3_df$consistency_rate_iwaba,
    level2_dominance_flips = sum(pair_df$flipped),
    stringsAsFactors = FALSE
  )
}

build_bf_df <- function(X, g, variables, alpha) {
  bf_df <- iwaba_brown_forsythe_matrix(X, g, varnames = variables, alpha = alpha)
  bf_df$analysis_family <- "brown_forsythe"
  bf_df$target_label <- bf_df$variable
  bf_df$estimate <- bf_df$statistic
  bf_df$ci_lower <- NA_real_
  bf_df$ci_upper <- NA_real_
  bf_df$ci_defined <- NA
  bf_df$boot_defined_rate <- NA_real_
  bf_df$reject_null <- !is.na(bf_df$p_value) & bf_df$p_value < alpha
  bf_df$notes <- "H0: equal within-group variances"
  bf_df
}

build_boot_ci_df <- function(X,
                             g,
                             variables,
                             radius_fn,
                             radius_name,
                             B,
                             seed) {
  pair_labels <- iwaba_pair_labels(combn(seq_along(variables), 2), variables)

  bh_boot_stat_fn <- function(Xb, gb) {
    result_b <- iwaba_full(Xb, gb, radius_fn = radius_fn, radius_name = radius_name)
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

  boot <- iwaba_within_group_bootstrap_vector(
    X = X,
    g = g,
    statistic_fn = bh_boot_stat_fn,
    B = B,
    conf_level = 0.95,
    seed = seed,
    min_finite = max(50L, floor(B * 0.50)),
    null_value = 0,
    progress = FALSE
  )

  boot_ci_df <- boot$ci
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
          "Bootstrap CI for signed heterogeneity-based radius correlation"
        )
      )
    )
  )

  list(
    ci = boot_ci_df,
    observed = boot$observed,
    B = boot$B
  )
}

build_branch_export <- function(branch_label,
                                scale_definition,
                                radius_definition,
                                entity_df,
                                pair_df,
                                level3_df,
                                level4_df,
                                level5_df,
                                summary_df,
                                bf_df,
                                boot_ci_df = NULL) {
  entity_export <- transform(
    entity_df,
    branch = branch_label,
    scale_definition = scale_definition,
    radius_definition = radius_definition,
    analysis_family = "entity",
    target_label = variable
  )

  pair_export <- transform(
    pair_df,
    branch = branch_label,
    scale_definition = scale_definition,
    radius_definition = radius_definition,
    analysis_family = "pairwise",
    target_label = pair_label
  )

  level3_export <- transform(
    level3_df,
    branch = branch_label,
    scale_definition = scale_definition,
    radius_definition = radius_definition,
    analysis_family = "level3",
    target_label = "global"
  )

  level4_export <- transform(
    level4_df,
    branch = branch_label,
    scale_definition = scale_definition,
    radius_definition = radius_definition,
    analysis_family = "level4",
    target_label = variable
  )

  level5_export <- transform(
    level5_df,
    branch = branch_label,
    scale_definition = scale_definition,
    radius_definition = radius_definition,
    analysis_family = "level5",
    target_label = pair_label
  )

  bf_export <- transform(
    bf_df,
    branch = branch_label,
    scale_definition = scale_definition,
    radius_definition = radius_definition
  )

  summary_export <- transform(
    summary_df,
    scale_definition = scale_definition,
    radius_definition = radius_definition,
    analysis_family = "summary",
    target_label = branch
  )

  boot_export <- NULL
  if (!is.null(boot_ci_df) && nrow(boot_ci_df) > 0L) {
    boot_export <- transform(
      boot_ci_df,
      branch = branch_label,
      scale_definition = scale_definition,
      radius_definition = radius_definition
    )
  }

  bind_named_rows(list(
    entity_export,
    pair_export,
    level3_export,
    level4_export,
    level5_export,
    bf_export,
    boot_export,
    summary_export
  ))
}

run_bh1996_branch <- function(X,
                              g,
                              variables,
                              branch_label,
                              scale_definition,
                              radius_fn,
                              radius_name,
                              radius_definition,
                              alpha,
                              B = 0L,
                              bootstrap_seed = NA_integer_) {
  result <- iwaba_full(X, g, radius_fn = radius_fn, radius_name = radius_name)
  entity_df <- build_entity_df(result, X, g, variables)
  pair_df <- build_pair_df(result, variables)
  level3_df <- build_level3_df(result, pair_df)
  level4_df <- build_level4_df(result, variables)
  level5_df <- build_level5_df(result, variables)
  summary_df <- build_summary_df(result, entity_df, pair_df, level3_df, level4_df, branch_label)
  bf_df <- build_bf_df(X, g, variables, alpha)

  boot_info <- NULL
  if (B > 0L && !is.na(bootstrap_seed)) {
    boot_info <- build_boot_ci_df(
      X = X,
      g = g,
      variables = variables,
      radius_fn = radius_fn,
      radius_name = radius_name,
      B = B,
      seed = bootstrap_seed
    )
  }

  export_df <- build_branch_export(
    branch_label = branch_label,
    scale_definition = scale_definition,
    radius_definition = radius_definition,
    entity_df = entity_df,
    pair_df = pair_df,
    level3_df = level3_df,
    level4_df = level4_df,
    level5_df = level5_df,
    summary_df = summary_df,
    bf_df = bf_df,
    boot_ci_df = if (is.null(boot_info)) NULL else boot_info$ci
  )

  list(
    result = result,
    entity = entity_df,
    pairwise = pair_df,
    level3 = level3_df,
    level4 = level4_df,
    level5 = level5_df,
    summary = summary_df,
    brown_forsythe = bf_df,
    bootstrap = boot_info,
    export = export_df,
    branch_label = branch_label,
    scale_definition = scale_definition,
    radius_definition = radius_definition
  )
}

write_primary_figure <- function(primary_branch, figure_filename, output_dir) {
  eta <- primary_branch$result$level1$eta_results
  l4 <- primary_branch$result$level4
  l5 <- primary_branch$result$level5
  pair_labels <- primary_branch$pairwise$pair_label
  pdf_path <- iwaba_open_pdf_device(figure_filename, output_dir = output_dir, width = 14, height = 10)

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

  boxplot(primary_branch$result$gq$radii,
          main = "(b) Within-Group SD Distributions",
          names = eta$variable, col = "#55A868",
          ylab = "Within-group SD",
          border = "darkgray", pch = 16, cex = 0.8)
  grid(col = "gray80", lty = 1)

  barplot(rbind(l4$E_C_int, l4$E_R_het_int), beside = TRUE,
          names.arg = eta$variable,
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
  pdf_path
}

write_radius_sensitivity_figure <- function(primary_branch,
                                            iqr_branch,
                                            figure_filename,
                                            output_dir) {
  pdf_path <- iwaba_open_pdf_device(figure_filename, output_dir = output_dir, width = 14, height = 10)

  par(mfrow = c(2, 2), mar = c(6, 5, 4, 2))

  eta_mat <- rbind(primary_branch$entity$eta_B_int, iqr_branch$entity$eta_B_int)
  colnames(eta_mat) <- primary_branch$entity$variable
  barplot(eta_mat, beside = TRUE,
          col = c("steelblue", "coral"),
          main = "(a) eta_B^(Int): SD vs IQR Radius",
          ylab = expression(eta[B]^"(Int)"),
          border = "darkgray",
          ylim = c(0, max(eta_mat) * 1.15))
  legend("topleft", legend = c("SD radius", "IQR radius"), fill = c("steelblue", "coral"), bty = "n")
  grid(col = "gray80", lty = 1)

  pi_mat <- rbind(primary_branch$level4$pi_R, iqr_branch$level4$pi_R)
  colnames(pi_mat) <- primary_branch$level4$variable
  barplot(pi_mat, beside = TRUE,
          col = c("steelblue", "coral"),
          main = expression("(b) " * pi[R] * ": SD vs IQR Radius"),
          ylab = expression(pi[R]),
          border = "darkgray",
          ylim = c(0, max(pi_mat) * 1.15))
  legend("topleft", legend = c("SD radius", "IQR radius"), fill = c("steelblue", "coral"), bty = "n")
  grid(col = "gray80", lty = 1)

  erhet_mat <- rbind(primary_branch$level4$E_R_het_int, iqr_branch$level4$E_R_het_int)
  colnames(erhet_mat) <- primary_branch$level4$variable
  barplot(erhet_mat, beside = TRUE,
          col = c("steelblue", "coral"),
          main = expression("(c) " * E[R*","*het]^"(Int): SD vs IQR Radius"),
          ylab = expression(E[R*","*het]^"(Int)"),
          border = "darkgray",
          ylim = c(0, max(erhet_mat) * 1.15))
  legend("topleft", legend = c("SD radius", "IQR radius"), fill = c("steelblue", "coral"), bty = "n")
  grid(col = "gray80", lty = 1)

  rhet_mat <- rbind(primary_branch$level5$r_R_het, iqr_branch$level5$r_R_het)
  colnames(rhet_mat) <- primary_branch$level5$pair_label
  y_min <- min(0, rhet_mat, na.rm = TRUE) * 1.15
  y_max <- max(rhet_mat, na.rm = TRUE) * 1.15
  barplot(rhet_mat, beside = TRUE,
          col = c("steelblue", "coral"),
          main = expression("(d) " * r[R*","*het] * ": SD vs IQR Radius"),
          ylab = expression(r[R*","*het]),
          border = "darkgray",
          ylim = c(y_min, y_max),
          names.arg = gsub("--", "\n", primary_branch$level5$pair_label, fixed = TRUE))
  abline(h = 0, col = "gray40", lwd = 1)
  legend("topleft", legend = c("SD radius", "IQR radius"), fill = c("steelblue", "coral"), bty = "n")
  grid(col = "gray80", lty = 1)

  dev.off()
  pdf_path
}

args <- commandArgs(trailingOnly = TRUE)
BOOT_B <- as.integer(parse_arg_value(args, "--B=", 399L))
ALPHA <- as.numeric(parse_arg_value(args, "--alpha=", 0.05))
BOOT_SEED <- as.integer(parse_arg_value(args, "--seed=", 20260321L))
OUTPUT_DIR <- parse_arg_value(args, "--output-dir=", "results")

PRIMARY_BOOT_SEED <- BOOT_SEED
IQR_BOOT_SEED <- BOOT_SEED + 1L

PRIMARY_FILE <- "bh1996_primary_v3.csv"
METADATA_FILE <- "bh1996_metadata_v1.csv"
PRIMARY_FIG_FILE <- "fig_bh1996_analysis_v3.pdf"
GROUP_SIZE_FILE <- "bh1996_group_sizes_v1.csv"
RADIUS_SENSITIVITY_FILE <- "bh1996_radius_sensitivity_v1.csv"
RADIUS_SENSITIVITY_FIG_FILE <- "fig_bh1996_radius_sensitivity_v1.pdf"
SCALE_SENSITIVITY_FILE <- "bh1996_scale_sensitivity_v1.csv"
INTERPRETIVE_FILE <- "bh1996_inference_interpretive_v1.csv"
RUNLOG_FILE <- file.path(OUTPUT_DIR, "realdata_runlog_v2.md")

SCRIPT_INVOCATION <- paste(c("real_data_bh1996.R", args), collapse = " ")
if (!nzchar(SCRIPT_INVOCATION)) {
  SCRIPT_INVOCATION <- "real_data_bh1996.R"
}

output_registry <- data.frame(
  label = character(),
  type = character(),
  path = character(),
  timestamp = character(),
  stringsAsFactors = FALSE
)

analysis_state <- list(
  started_at = timestamp_now(),
  primary_summary = NULL,
  iqr_summary = NULL,
  standardized_summary = NULL,
  group_size_summary = NULL
)

update_metadata <- function() {
  metadata_list <- list(
    script = "real_data_bh1996.R",
    dataset = "bh1996",
    variables = "COHES; LEAD; WBEING; HRS",
    grouping_variable = "GRP",
    primary_radius_definition = "Population within-company SD",
    sensitivity_radius_definition = "0.7413 * IQR within company",
    scale_sensitivity_definition = "Global z-score by variable over all individuals",
    bootstrap_resamples_primary = BOOT_B,
    bootstrap_resamples_iqr = BOOT_B,
    alpha = ALPHA,
    primary_boot_seed = PRIMARY_BOOT_SEED,
    iqr_boot_seed = IQR_BOOT_SEED,
    started_at = analysis_state$started_at,
    metadata_updated_at = timestamp_now(),
    output_timestamps = format_output_registry(output_registry),
    outputs_created = if (nrow(output_registry) == 0L) "" else paste(output_registry$path, collapse = "; ")
  )

  if (!is.null(analysis_state$group_size_summary)) {
    metadata_list$K <- analysis_state$group_size_summary$K
    metadata_list$n <- analysis_state$group_size_summary$n
    metadata_list$min_n_k <- analysis_state$group_size_summary$min_n_k
    metadata_list$q1_n_k <- analysis_state$group_size_summary$q1_n_k
    metadata_list$median_n_k <- analysis_state$group_size_summary$median_n_k
    metadata_list$mean_n_k <- analysis_state$group_size_summary$mean_n_k
    metadata_list$q3_n_k <- analysis_state$group_size_summary$q3_n_k
    metadata_list$max_n_k <- analysis_state$group_size_summary$max_n_k
  }

  if (!is.null(analysis_state$primary_summary)) {
    metadata_list$primary_G <- analysis_state$primary_summary$G
    metadata_list$primary_radius_prop <- analysis_state$primary_summary$radius_prop
    metadata_list$primary_mean_E_R_het_int <- analysis_state$primary_summary$mean_E_R_het_int
  }

  if (!is.null(analysis_state$iqr_summary)) {
    metadata_list$iqr_G <- analysis_state$iqr_summary$G
    metadata_list$iqr_radius_prop <- analysis_state$iqr_summary$radius_prop
    metadata_list$iqr_mean_E_R_het_int <- analysis_state$iqr_summary$mean_E_R_het_int
  }

  if (!is.null(analysis_state$standardized_summary)) {
    metadata_list$standardized_G <- analysis_state$standardized_summary$G
    metadata_list$standardized_radius_prop <- analysis_state$standardized_summary$radius_prop
    metadata_list$standardized_mean_E_R_het_int <- analysis_state$standardized_summary$mean_E_R_het_int
  }

  metadata_paths <- iwaba_save_metadata(
    metadata = metadata_list,
    metadata_filename = METADATA_FILE,
    output_dir = OUTPUT_DIR
  )

  output_registry <<- output_registry[output_registry$label != "metadata", , drop = FALSE]
  output_registry <<- register_dual_output(
    output_registry,
    "metadata",
    list(csv = metadata_paths$csv, rds = metadata_paths$rds, timestamp = timestamp_now())
  )
  metadata_paths
}

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
X_standardized <- scale(X)
X_standardized <- as.matrix(X_standardized)
colnames(X_standardized) <- variables

cat(sprintf("  n = %d, K = %d, p = %d\n", nrow(X), nlevels(g), ncol(X)))
cat(sprintf("  Primary bootstrap B = %d, alpha = %.3f, primary seed = %d, IQR seed = %d\n",
            BOOT_B, ALPHA, PRIMARY_BOOT_SEED, IQR_BOOT_SEED))

cat("\nSubtask A: Primary SD-radius rerun\n")
primary_branch <- run_bh1996_branch(
  X = X,
  g = g,
  variables = variables,
  branch_label = "original_scale_sd_radius",
  scale_definition = "Original scale",
  radius_fn = iwaba_radius_sd,
  radius_name = "sd",
  radius_definition = "Population within-company SD",
  alpha = ALPHA,
  B = BOOT_B,
  bootstrap_seed = PRIMARY_BOOT_SEED
)
analysis_state$primary_summary <- primary_branch$summary

primary_paths <- save_dual_output(
  csv_object = primary_branch$export,
  rds_object = primary_branch,
  csv_filename = PRIMARY_FILE,
  output_dir = OUTPUT_DIR
)
output_registry <- register_dual_output(output_registry, "primary", primary_paths)

primary_fig_path <- write_primary_figure(primary_branch, PRIMARY_FIG_FILE, OUTPUT_DIR)
output_registry <- register_pdf_output(output_registry, "primary_figure", primary_fig_path)

metadata_paths <- update_metadata()

append_runlog_entry(
  runlog_file = RUNLOG_FILE,
  subtask = "Task A - Primary SD-radius rerun",
  files_changed = c(basename(primary_paths$csv), basename(primary_paths$rds), basename(primary_fig_path), basename(metadata_paths$csv), basename(metadata_paths$rds)),
  commands_run = c(
    SCRIPT_INVOCATION,
    sprintf("Subtask A: original-scale SD-radius analysis with within-group bootstrap B=%d", BOOT_B)
  ),
  outputs_created = c(primary_paths$csv, primary_paths$rds, primary_fig_path, metadata_paths$csv, metadata_paths$rds),
  seeds_used = sprintf("Primary bootstrap seed = %d", PRIMARY_BOOT_SEED),
  what_remains = c(
    "Task B: group-size summary",
    "Task C: IQR-radius sensitivity",
    "Task D: scale-comparability sensitivity",
    "Task E: interpretive inference table"
  )
)

cat("\nSubtask B: Group-size summary\n")
group_size_summary <- build_group_size_summary(g)
analysis_state$group_size_summary <- group_size_summary

group_size_paths <- save_dual_output(
  csv_object = group_size_summary,
  rds_object = list(summary = group_size_summary, company_sizes = as.numeric(table(g))),
  csv_filename = GROUP_SIZE_FILE,
  output_dir = OUTPUT_DIR
)
output_registry <- register_dual_output(output_registry, "group_sizes", group_size_paths)

metadata_paths <- update_metadata()

append_runlog_entry(
  runlog_file = RUNLOG_FILE,
  subtask = "Task B - Group-size summary",
  files_changed = c(basename(group_size_paths$csv), basename(group_size_paths$rds), basename(metadata_paths$csv), basename(metadata_paths$rds)),
  commands_run = c(
    SCRIPT_INVOCATION,
    "Subtask B: company-size summary from table(GRP)"
  ),
  outputs_created = c(group_size_paths$csv, group_size_paths$rds, metadata_paths$csv, metadata_paths$rds),
  seeds_used = "No RNG used",
  what_remains = c(
    "Task C: IQR-radius sensitivity",
    "Task D: scale-comparability sensitivity",
    "Task E: interpretive inference table"
  )
)

cat("\nSubtask C: IQR-radius sensitivity\n")
iqr_branch <- run_bh1996_branch(
  X = X,
  g = g,
  variables = variables,
  branch_label = "original_scale_iqr_radius",
  scale_definition = "Original scale",
  radius_fn = iwaba_radius_iqr,
  radius_name = "iqr",
  radius_definition = "0.7413 * IQR within company",
  alpha = ALPHA,
  B = BOOT_B,
  bootstrap_seed = IQR_BOOT_SEED
)
analysis_state$iqr_summary <- iqr_branch$summary

radius_sensitivity_paths <- save_dual_output(
  csv_object = iqr_branch$export,
  rds_object = iqr_branch,
  csv_filename = RADIUS_SENSITIVITY_FILE,
  output_dir = OUTPUT_DIR
)
output_registry <- register_dual_output(output_registry, "radius_sensitivity", radius_sensitivity_paths)

radius_sensitivity_fig_path <- write_radius_sensitivity_figure(
  primary_branch = primary_branch,
  iqr_branch = iqr_branch,
  figure_filename = RADIUS_SENSITIVITY_FIG_FILE,
  output_dir = OUTPUT_DIR
)
output_registry <- register_pdf_output(output_registry, "radius_sensitivity_figure", radius_sensitivity_fig_path)

metadata_paths <- update_metadata()

append_runlog_entry(
  runlog_file = RUNLOG_FILE,
  subtask = "Task C - IQR-radius sensitivity",
  files_changed = c(basename(radius_sensitivity_paths$csv), basename(radius_sensitivity_paths$rds), basename(radius_sensitivity_fig_path), basename(metadata_paths$csv), basename(metadata_paths$rds)),
  commands_run = c(
    SCRIPT_INVOCATION,
    sprintf("Subtask C: original-scale IQR-radius analysis with within-group bootstrap B=%d", BOOT_B)
  ),
  outputs_created = c(radius_sensitivity_paths$csv, radius_sensitivity_paths$rds, radius_sensitivity_fig_path, metadata_paths$csv, metadata_paths$rds),
  seeds_used = sprintf("IQR bootstrap seed = %d", IQR_BOOT_SEED),
  what_remains = c(
    "Task D: scale-comparability sensitivity",
    "Task E: interpretive inference table"
  )
)

cat("\nSubtask D: Scale-comparability sensitivity\n")
standardized_branch <- run_bh1996_branch(
  X = X_standardized,
  g = g,
  variables = variables,
  branch_label = "standardized_scale_sd_radius",
  scale_definition = "Global z-score by variable",
  radius_fn = iwaba_radius_sd,
  radius_name = "sd",
  radius_definition = "Population within-company SD on globally standardized variables",
  alpha = ALPHA,
  B = 0L,
  bootstrap_seed = NA_integer_
)
analysis_state$standardized_summary <- standardized_branch$summary

scale_sensitivity_df <- rbind(
  transform(
    merge(primary_branch$entity[, c("variable", "eta_B_waba", "eta_B_int", "E_waba", "E_int", "V_C", "V_R", "V_W")],
          primary_branch$level4[, c("variable", "pi_R", "E_R_het_int")],
          by = "variable"),
    branch = "original_scale_sd_radius",
    scale_definition = "Original scale",
    radius_definition = "Population within-company SD"
  ),
  transform(
    merge(standardized_branch$entity[, c("variable", "eta_B_waba", "eta_B_int", "E_waba", "E_int", "V_C", "V_R", "V_W")],
          standardized_branch$level4[, c("variable", "pi_R", "E_R_het_int")],
          by = "variable"),
    branch = "standardized_scale_sd_radius",
    scale_definition = "Global z-score by variable",
    radius_definition = "Population within-company SD"
  ),
  transform(
    merge(iqr_branch$entity[, c("variable", "eta_B_waba", "eta_B_int", "E_waba", "E_int", "V_C", "V_R", "V_W")],
          iqr_branch$level4[, c("variable", "pi_R", "E_R_het_int")],
          by = "variable"),
    branch = "original_scale_iqr_radius",
    scale_definition = "Original scale",
    radius_definition = "0.7413 * IQR within company"
  )
)

scale_sensitivity_paths <- save_dual_output(
  csv_object = scale_sensitivity_df,
  rds_object = list(
    comparison = scale_sensitivity_df,
    standardized_branch = standardized_branch
  ),
  csv_filename = SCALE_SENSITIVITY_FILE,
  output_dir = OUTPUT_DIR
)
output_registry <- register_dual_output(output_registry, "scale_sensitivity", scale_sensitivity_paths)

metadata_paths <- update_metadata()

append_runlog_entry(
  runlog_file = RUNLOG_FILE,
  subtask = "Task D - Scale-comparability sensitivity",
  files_changed = c(basename(scale_sensitivity_paths$csv), basename(scale_sensitivity_paths$rds), basename(metadata_paths$csv), basename(metadata_paths$rds)),
  commands_run = c(
    SCRIPT_INVOCATION,
    "Subtask D: globally standardized SD-radius rerun without bootstrap"
  ),
  outputs_created = c(scale_sensitivity_paths$csv, scale_sensitivity_paths$rds, metadata_paths$csv, metadata_paths$rds),
  seeds_used = "No RNG used",
  what_remains = "Task E: interpretive inference table"
)

cat("\nSubtask E: Interpretive inference table\n")
primary_boot_delta <- primary_branch$bootstrap$ci
primary_boot_delta <- primary_boot_delta[primary_boot_delta$analysis_family == "delta_E", , drop = FALSE]
primary_boot_delta <- primary_boot_delta[, c("target_label", "ci_lower", "ci_upper", "ci_defined", "boot_defined_rate"), drop = FALSE]
names(primary_boot_delta)[names(primary_boot_delta) == "target_label"] <- "variable"

interpretive_df <- merge(
  merge(
    primary_branch$brown_forsythe[, c("variable", "p_value", "p_adj", "reject"), drop = FALSE],
    primary_branch$level4[, c("variable", "pi_R", "E_R_het_int"), drop = FALSE],
    by = "variable"
  ),
  primary_branch$entity[, c("variable", "class_iwaba", "delta_E"), drop = FALSE],
  by = "variable"
)
interpretive_df <- merge(interpretive_df, primary_boot_delta, by = "variable", all.x = TRUE)
interpretive_df <- interpretive_df[match(variables, interpretive_df$variable), , drop = FALSE]
interpretive_df$within_group_dominant <- interpretive_df$class_iwaba == "Within"
names(interpretive_df)[names(interpretive_df) == "p_value"] <- "bf_p_value"
names(interpretive_df)[names(interpretive_df) == "p_adj"] <- "bf_p_adj"
names(interpretive_df)[names(interpretive_df) == "reject"] <- "bf_bh_reject"

interpretive_paths <- save_dual_output(
  csv_object = interpretive_df,
  rds_object = interpretive_df,
  csv_filename = INTERPRETIVE_FILE,
  output_dir = OUTPUT_DIR
)
output_registry <- register_dual_output(output_registry, "interpretive", interpretive_paths)

metadata_paths <- update_metadata()

append_runlog_entry(
  runlog_file = RUNLOG_FILE,
  subtask = "Task E - Interpretive inference table",
  files_changed = c(basename(interpretive_paths$csv), basename(interpretive_paths$rds), basename(metadata_paths$csv), basename(metadata_paths$rds)),
  commands_run = c(
    SCRIPT_INVOCATION,
    "Subtask E: merge Brown-Forsythe, Level IV, and primary bootstrap delta_E intervals"
  ),
  outputs_created = c(interpretive_paths$csv, interpretive_paths$rds, metadata_paths$csv, metadata_paths$rds),
  seeds_used = sprintf("Primary bootstrap seed reused for interpretive delta_E intervals = %d", PRIMARY_BOOT_SEED),
  what_remains = "bh1996 milestone complete"
)

cat("\nSaved outputs:\n")
for (path in output_registry$path) {
  cat(sprintf("  %s\n", path))
}

cat("\nKey summaries:\n")
cat(sprintf("  Primary: G = %.3f, radius share = %.3f, mean E_R,het = %.3f\n",
            analysis_state$primary_summary$G,
            analysis_state$primary_summary$radius_prop,
            analysis_state$primary_summary$mean_E_R_het_int))
cat(sprintf("  IQR sensitivity: G = %.3f, radius share = %.3f, mean E_R,het = %.3f\n",
            analysis_state$iqr_summary$G,
            analysis_state$iqr_summary$radius_prop,
            analysis_state$iqr_summary$mean_E_R_het_int))
cat(sprintf("  Standardized: G = %.3f, radius share = %.3f, mean E_R,het = %.3f\n",
            analysis_state$standardized_summary$G,
            analysis_state$standardized_summary$radius_prop,
            analysis_state$standardized_summary$mean_E_R_het_int))

cat("\nDone.\n")
