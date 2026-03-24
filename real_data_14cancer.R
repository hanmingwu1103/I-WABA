################################################################################
## I-WABA Real Data Analysis: 14-Cancer Gene Expression Dataset
## Screening and selection-sensitivity milestone
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
    sprintf("## Milestone: 14-Cancer / %s", subtask),
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
                             tag = NULL,
                             row.names = FALSE) {
  csv_path <- iwaba_output_path(csv_filename, output_dir = output_dir, tag = tag, create_dir = TRUE)
  rds_path <- iwaba_output_path(
    sub("\\.csv$", ".rds", csv_filename, ignore.case = TRUE),
    output_dir = output_dir,
    tag = tag,
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

register_single_output <- function(output_registry, label, path, type = "file") {
  rbind(
    output_registry,
    data.frame(
      label = label,
      type = type,
      path = path,
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

build_selected_genes_df <- function(rank_df, rank_positions, branch_label) {
  out <- rank_df[rank_positions, c("variable", "F_score", "rank"), drop = FALSE]
  out$branch <- branch_label
  out[order(out$rank), c("branch", "variable", "F_score", "rank"), drop = FALSE]
}

build_level1_df <- function(result_or_level1, varnames, f_scores_lookup) {
  l1 <- if (!is.null(result_or_level1$level1)) result_or_level1$level1 else result_or_level1
  eta <- l1$eta_results
  delta_E <- iwaba_entity_dominance_contrast(result_or_level1)
  names(delta_E) <- varnames

  data.frame(
    variable = varnames,
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
    delta_E = unname(delta_E[varnames]),
    F_score = unname(f_scores_lookup[varnames]),
    stringsAsFactors = FALSE
  )
}

build_level2_df <- function(result, varnames) {
  l2 <- result$level2
  pair_labels <- iwaba_pair_labels(l2$pairs, varnames)
  delta_A <- iwaba_pairwise_dominance_contrast(result, varnames = varnames)

  data.frame(
    pair_label = pair_labels,
    var1 = varnames[l2$pairs[1, ]],
    var2 = varnames[l2$pairs[2, ]],
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

build_level2_aggregate_df <- function(result, level2_df) {
  data.frame(
    n_pairs = nrow(level2_df),
    mean_abs_r_between_waba = mean(abs(level2_df$r_between_waba)),
    mean_abs_r_between_int = mean(abs(level2_df$r_between_int)),
    mean_abs_weighted_between_waba = mean(abs(level2_df$weighted_between_waba)),
    mean_abs_weighted_within_waba = mean(abs(level2_df$weighted_within_waba)),
    mean_abs_weighted_between_int = mean(abs(level2_df$weighted_between_int)),
    mean_abs_weighted_within_int = mean(abs(level2_df$weighted_within_int)),
    proportion_between_dominant_waba = mean(level2_df$dominance_waba == "Between"),
    proportion_between_dominant_iwaba = mean(level2_df$dominance_iwaba == "Between"),
    consistency_rate_waba = result$level3$consistency_rate_waba,
    consistency_rate_iwaba = result$level3$consistency_rate_int,
    stringsAsFactors = FALSE
  )
}

build_level3_df <- function(result, level2_df) {
  data.frame(
    consistency_rate_waba = result$level3$consistency_rate_waba,
    consistency_rate_iwaba = result$level3$consistency_rate_int,
    n_pairs = nrow(level2_df),
    n_dominance_flips = sum(level2_df$flipped),
    stringsAsFactors = FALSE
  )
}

build_level4_df <- function(result_or_level4, varnames) {
  l4 <- if (!is.null(result_or_level4$level4)) result_or_level4$level4 else result_or_level4

  data.frame(
    variable = varnames,
    pi_R = l4$pi_R,
    E_C_int = l4$E_C_int,
    E_R_int = l4$E_R_int,
    E_R_het_int = l4$E_R_het_int,
    V_R_het = l4$V_R_het,
    stringsAsFactors = FALSE
  )
}

build_level5_df <- function(result, varnames) {
  l2 <- result$level2
  l5 <- result$level5
  pair_labels <- iwaba_pair_labels(l2$pairs, varnames)

  data.frame(
    pair_label = pair_labels,
    var1 = varnames[l2$pairs[1, ]],
    var2 = varnames[l2$pairs[2, ]],
    r_R_billard = l5$r_R_billard,
    r_R_billard_defined = l5$r_R_billard_defined,
    r_R_het = l5$r_R_het,
    r_R_het_defined = l5$r_R_het_defined,
    stringsAsFactors = FALSE
  )
}

build_screening_df <- function(X, g, entity_df, level4_df, alpha) {
  screening_df <- iwaba_brown_forsythe_matrix(
    X,
    g,
    varnames = entity_df$variable,
    alpha = alpha
  )

  screening_df$F_score <- entity_df$F_score[match(screening_df$variable, entity_df$variable)]
  screening_df$eta_B_int <- entity_df$eta_B_int[match(screening_df$variable, entity_df$variable)]
  screening_df$delta_eta_B <- entity_df$delta_eta_B[match(screening_df$variable, entity_df$variable)]
  screening_df$pi_R <- level4_df$pi_R[match(screening_df$variable, level4_df$variable)]
  screening_df$E_R_het_int <- level4_df$E_R_het_int[match(screening_df$variable, level4_df$variable)]
  screening_df$within_group_dominant <- entity_df$class_iwaba[match(screening_df$variable, entity_df$variable)] == "Within"
  screening_df
}
build_selected_summary_df <- function(result,
                                      pair_df,
                                      pair_agg_df,
                                      screening_df,
                                      branch_label,
                                      selection_rule,
                                      bootstrap_B,
                                      bootstrap_seed,
                                      selected_pair_count) {
  l1 <- result$level1
  l4 <- result$level4
  l5 <- result$level5

  data.frame(
    branch = branch_label,
    selection_rule = selection_rule,
    n = l1$n,
    p = l1$p,
    K = l1$K,
    G = l1$info_gain,
    radius_prop = l1$radius_prop,
    mean_eta_B_waba = mean(l1$eta_results$eta_B_waba),
    mean_eta_B_int = mean(l1$eta_results$eta_B_iwaba),
    mean_delta_eta_B = mean(l1$eta_results$delta_eta),
    reclassified_count = sum(l1$eta_results$class_waba != l1$eta_results$class_iwaba),
    reclassified_prop = mean(l1$eta_results$class_waba != l1$eta_results$class_iwaba),
    mean_pi_R = mean(l4$pi_R),
    mean_E_R_het_int = mean(l4$E_R_het_int),
    T_het = iwaba_global_t_het(result),
    mean_abs_r_between_waba = pair_agg_df$mean_abs_r_between_waba,
    mean_abs_r_between_int = pair_agg_df$mean_abs_r_between_int,
    mean_abs_weighted_between_waba = pair_agg_df$mean_abs_weighted_between_waba,
    mean_abs_weighted_between_int = pair_agg_df$mean_abs_weighted_between_int,
    mean_abs_weighted_within_waba = pair_agg_df$mean_abs_weighted_within_waba,
    mean_abs_weighted_within_int = pair_agg_df$mean_abs_weighted_within_int,
    proportion_between_dominant_waba = pair_agg_df$proportion_between_dominant_waba,
    proportion_between_dominant_iwaba = pair_agg_df$proportion_between_dominant_iwaba,
    consistency_rate_waba = result$level3$consistency_rate_waba,
    consistency_rate_iwaba = result$level3$consistency_rate_int,
    dominance_flip_count = sum(pair_df$flipped),
    mean_r_R_billard = l5$mean_r_R_billard,
    mean_r_R_het_signed = l5$mean_r_R_het_signed,
    positive_r_R_het_count = sum(l5$r_R_het > 0, na.rm = TRUE),
    negative_r_R_het_count = sum(l5$r_R_het < 0, na.rm = TRUE),
    bf_raw_count = sum(screening_df$p_value < ALPHA, na.rm = TRUE),
    bf_bh_count = sum(screening_df$reject, na.rm = TRUE),
    bf_bh_prop = mean(screening_df$reject, na.rm = TRUE),
    selected_pair_count = selected_pair_count,
    bootstrap_B = bootstrap_B,
    bootstrap_seed = bootstrap_seed,
    stringsAsFactors = FALSE
  )
}

build_allgenes_summary_df <- function(level1_obj,
                                      level4_obj,
                                      entity_df,
                                      screening_df,
                                      selection_rule) {
  mean_radius_all <- rowMeans(level1_obj$gq$radii)

  data.frame(
    branch = "allgenes",
    selection_rule = selection_rule,
    n = level1_obj$n,
    p = level1_obj$p,
    K = level1_obj$K,
    G = level1_obj$info_gain,
    radius_prop = level1_obj$radius_prop,
    mean_eta_B_waba = mean(entity_df$eta_B_waba),
    mean_eta_B_int = mean(entity_df$eta_B_int),
    mean_delta_eta_B = mean(entity_df$delta_eta_B),
    reclassified_count = sum(entity_df$class_waba != entity_df$class_iwaba),
    reclassified_prop = mean(entity_df$class_waba != entity_df$class_iwaba),
    mean_pi_R = mean(level4_obj$pi_R),
    mean_E_R_het_int = mean(level4_obj$E_R_het_int),
    T_het = iwaba_global_t_het(level4_obj),
    mean_abs_r_between_waba = NA_real_,
    mean_abs_r_between_int = NA_real_,
    mean_abs_weighted_between_waba = NA_real_,
    mean_abs_weighted_between_int = NA_real_,
    mean_abs_weighted_within_waba = NA_real_,
    mean_abs_weighted_within_int = NA_real_,
    proportion_between_dominant_waba = NA_real_,
    proportion_between_dominant_iwaba = NA_real_,
    consistency_rate_waba = NA_real_,
    consistency_rate_iwaba = NA_real_,
    dominance_flip_count = NA_real_,
    mean_r_R_billard = NA_real_,
    mean_r_R_het_signed = NA_real_,
    positive_r_R_het_count = NA_real_,
    negative_r_R_het_count = NA_real_,
    bf_raw_count = sum(screening_df$p_value < ALPHA, na.rm = TRUE),
    bf_bh_count = sum(screening_df$reject, na.rm = TRUE),
    bf_bh_prop = mean(screening_df$reject, na.rm = TRUE),
    selected_pair_count = 0,
    bootstrap_B = 0L,
    bootstrap_seed = NA_integer_,
    mean_within_group_sd_min = min(mean_radius_all),
    mean_within_group_sd_max = max(mean_radius_all),
    mean_within_group_sd_ratio = max(mean_radius_all) / min(mean_radius_all),
    stringsAsFactors = FALSE
  )
}

build_selected_branch_export <- function(branch_label,
                                         selection_rule,
                                         entity_df,
                                         pair_df,
                                         pair_agg_df,
                                         level3_df,
                                         level4_df,
                                         level5_df,
                                         screening_df,
                                         summary_df,
                                         boot_ci_df = NULL) {
  entity_export <- transform(entity_df, branch = branch_label, selection_rule = selection_rule,
                             analysis_family = "entity", target_label = variable)
  pair_export <- transform(pair_df, branch = branch_label, selection_rule = selection_rule,
                           analysis_family = "pairwise", target_label = pair_label)
  pair_agg_export <- transform(pair_agg_df, branch = branch_label, selection_rule = selection_rule,
                               analysis_family = "pairwise_aggregate", target_label = "all_pairs")
  level3_export <- transform(level3_df, branch = branch_label, selection_rule = selection_rule,
                             analysis_family = "level3", target_label = "global")
  level4_export <- transform(level4_df, branch = branch_label, selection_rule = selection_rule,
                             analysis_family = "level4", target_label = variable)
  level5_export <- transform(level5_df, branch = branch_label, selection_rule = selection_rule,
                             analysis_family = "level5", target_label = pair_label)
  screening_export <- transform(screening_df, branch = branch_label, selection_rule = selection_rule,
                                analysis_family = "brown_forsythe", target_label = variable)
  boot_export <- NULL
  if (!is.null(boot_ci_df) && nrow(boot_ci_df) > 0L) {
    boot_export <- transform(boot_ci_df, branch = branch_label, selection_rule = selection_rule)
  }
  summary_export <- transform(summary_df, branch = branch_label, selection_rule = selection_rule,
                              analysis_family = "summary", target_label = branch)

  bind_named_rows(list(
    entity_export,
    pair_export,
    pair_agg_export,
    level3_export,
    level4_export,
    level5_export,
    screening_export,
    boot_export,
    summary_export
  ))
}

build_allgenes_export <- function(entity_df,
                                  level4_df,
                                  screening_df,
                                  summary_df,
                                  selection_rule) {
  entity_export <- transform(entity_df, branch = "allgenes", selection_rule = selection_rule,
                             analysis_family = "entity", target_label = variable)
  level4_export <- transform(level4_df, branch = "allgenes", selection_rule = selection_rule,
                             analysis_family = "level4", target_label = variable)
  screening_export <- transform(screening_df, branch = "allgenes", selection_rule = selection_rule,
                                analysis_family = "brown_forsythe", target_label = variable)
  summary_export <- transform(summary_df, branch = "allgenes", selection_rule = selection_rule,
                              analysis_family = "summary", target_label = "allgenes")

  bind_named_rows(list(entity_export, level4_export, screening_export, summary_export))
}

run_selected_branch <- function(X,
                                g,
                                branch_label,
                                selection_rule,
                                f_scores_lookup,
                                alpha,
                                B,
                                bootstrap_seed,
                                max_pair_boot) {
  result <- iwaba_full(X, g)
  varnames <- colnames(X)

  entity_df <- build_level1_df(result, varnames, f_scores_lookup)
  pair_df <- build_level2_df(result, varnames)
  pair_agg_df <- build_level2_aggregate_df(result, pair_df)
  level3_df <- build_level3_df(result, pair_df)
  level4_df <- build_level4_df(result, varnames)
  level5_df <- build_level5_df(result, varnames)
  screening_df <- build_screening_df(X, g, entity_df, level4_df, alpha)

  screened_genes <- screening_df$variable[screening_df$reject]
  candidate_pair_idx <- if (length(screened_genes) >= 2L) {
    which(pair_df$var1 %in% screened_genes & pair_df$var2 %in% screened_genes)
  } else {
    integer(0)
  }
  if (length(candidate_pair_idx) == 0L) {
    candidate_pair_idx <- seq_len(nrow(pair_df))
  }
  candidate_pair_idx <- candidate_pair_idx[order(abs(level5_df$r_R_het[candidate_pair_idx]), decreasing = TRUE)]
  selected_pair_idx <- head(candidate_pair_idx, min(max_pair_boot, length(candidate_pair_idx)))
  selected_pair_labels <- pair_df$pair_label[selected_pair_idx]

  boot_info <- NULL
  if (B > 0L) {
    selected_pair_labels_local <- selected_pair_labels
    varnames_local <- varnames

    boot_stat_fn <- function(Xb, gb) {
      result_b <- iwaba_full(Xb, gb)
      pair_labels_b <- iwaba_pair_labels(result_b$level2$pairs, varnames_local)
      delta_A_b <- iwaba_pairwise_dominance_contrast(result_b, varnames = varnames_local)
      names(delta_A_b) <- pair_labels_b
      r_R_het_b <- result_b$level5$r_R_het
      names(r_R_het_b) <- pair_labels_b

      stats::setNames(
        c(
          mean(result_b$level1$eta_results$eta_B_iwaba),
          mean(result_b$level4$pi_R),
          mean(result_b$level4$E_R_het_int),
          iwaba_global_t_het(result_b),
          delta_A_b[selected_pair_labels_local],
          r_R_het_b[selected_pair_labels_local]
        ),
        c(
          "aggregate::mean_eta_B_int",
          "aggregate::mean_pi_R",
          "aggregate::mean_E_R_het_int",
          "aggregate::T_het",
          paste0("delta_A::", selected_pair_labels_local),
          paste0("r_R_het::", selected_pair_labels_local)
        )
      )
    }

    boot_info <- iwaba_within_group_bootstrap_vector(
      X = X,
      g = g,
      statistic_fn = boot_stat_fn,
      B = B,
      conf_level = 0.95,
      seed = bootstrap_seed,
      min_finite = max(50L, floor(B * 0.50)),
      null_value = 0,
      progress = FALSE
    )

    boot_ci_df <- boot_info$ci
    name_parts <- strsplit(boot_ci_df$stat_name, "::", fixed = TRUE)
    boot_ci_df$analysis_family <- vapply(name_parts, `[`, character(1), 1)
    boot_ci_df$target_label <- vapply(name_parts, `[`, character(1), 2)
    boot_ci_df$p_value <- NA_real_
    boot_ci_df$p_adj <- NA_real_
    boot_ci_df$reject <- NA
    boot_ci_df$notes <- ifelse(
      boot_ci_df$analysis_family == "aggregate",
      sprintf("Within-group bootstrap CI for aggregate %s summary", branch_label),
      ifelse(
        boot_ci_df$analysis_family == "delta_A",
        sprintf("Restricted pairwise bootstrap CI for %s pairwise dominance contrast", branch_label),
        sprintf("Restricted pairwise bootstrap CI for %s signed radius heterogeneity correlation", branch_label)
      )
    )
    boot_info$ci <- boot_ci_df
  }

  summary_df <- build_selected_summary_df(
    result = result,
    pair_df = pair_df,
    pair_agg_df = pair_agg_df,
    screening_df = screening_df,
    branch_label = branch_label,
    selection_rule = selection_rule,
    bootstrap_B = B,
    bootstrap_seed = bootstrap_seed,
    selected_pair_count = length(selected_pair_labels)
  )

  export_df <- build_selected_branch_export(
    branch_label = branch_label,
    selection_rule = selection_rule,
    entity_df = entity_df,
    pair_df = pair_df,
    pair_agg_df = pair_agg_df,
    level3_df = level3_df,
    level4_df = level4_df,
    level5_df = level5_df,
    screening_df = screening_df,
    summary_df = summary_df,
    boot_ci_df = if (is.null(boot_info)) NULL else boot_info$ci
  )

  list(
    result = result,
    entity = entity_df,
    pairwise = pair_df,
    pairwise_aggregate = pair_agg_df,
    level3 = level3_df,
    level4 = level4_df,
    level5 = level5_df,
    screening = screening_df,
    bootstrap = boot_info,
    summary = summary_df,
    export = export_df,
    selected_pair_labels = selected_pair_labels,
    selected_pair_table = pair_df[selected_pair_idx, , drop = FALSE],
    branch = branch_label,
    selection_rule = selection_rule
  )
}
run_allgenes_branch <- function(X, g, f_scores_lookup, alpha, selection_rule) {
  level1_obj <- iwaba_level1(X, g)
  level4_obj <- iwaba_level4(level1_obj$gq, level1_obj$V_C, level1_obj$V_R, level1_obj$V_W)

  entity_df <- build_level1_df(list(level1 = level1_obj), colnames(X), f_scores_lookup)
  level4_df <- build_level4_df(list(level4 = level4_obj), colnames(X))
  screening_df <- build_screening_df(X, g, entity_df, level4_df, alpha)
  summary_df <- build_allgenes_summary_df(
    level1_obj = level1_obj,
    level4_obj = level4_obj,
    entity_df = entity_df,
    screening_df = screening_df,
    selection_rule = selection_rule
  )
  export_df <- build_allgenes_export(
    entity_df = entity_df,
    level4_df = level4_df,
    screening_df = screening_df,
    summary_df = summary_df,
    selection_rule = selection_rule
  )

  list(
    level1 = level1_obj,
    level4 = level4_obj,
    entity = entity_df,
    screening = screening_df,
    summary = summary_df,
    export = export_df,
    branch = "allgenes",
    selection_rule = selection_rule
  )
}

bootstrap_two_sided_p <- function(samples, null_value = 0, min_finite = 50L) {
  finite_samples <- samples[is.finite(samples)]
  n_finite <- length(finite_samples)

  if (n_finite < min_finite) {
    return(NA_real_)
  }

  lower_tail <- (sum(finite_samples <= null_value) + 1) / (n_finite + 1)
  upper_tail <- (sum(finite_samples >= null_value) + 1) / (n_finite + 1)
  min(1, 2 * min(lower_tail, upper_tail))
}

checkpointed_within_group_bootstrap_vector <- function(X, g, statistic_fn,
                                                       B = 999L,
                                                       conf_level = 0.95,
                                                       seed = NULL,
                                                       min_finite = 50L,
                                                       null_value = NULL,
                                                       checkpoint_file = NULL,
                                                       checkpoint_every = 25L,
                                                       progress = FALSE) {
  X <- as.matrix(X)
  g <- as.factor(g)
  checkpoint_every <- max(1L, as.integer(checkpoint_every))

  checkpoint_exists <- !is.null(checkpoint_file) && file.exists(checkpoint_file)
  resumed_from_checkpoint <- FALSE

  if (checkpoint_exists) {
    checkpoint <- readRDS(checkpoint_file)
    if (is.null(checkpoint$B) || checkpoint$B != B) {
      stop("Checkpoint bootstrap size does not match requested B.")
    }

    observed <- checkpoint$observed
    stat_names <- checkpoint$stat_names
    boot_matrix <- checkpoint$boot_matrix
    completed <- checkpoint$completed
    resumed_from_checkpoint <- completed > 0L

    if (!is.null(checkpoint$rng_state)) {
      assign(".Random.seed", checkpoint$rng_state, envir = .GlobalEnv)
    } else if (!is.null(seed)) {
      set.seed(seed)
    }
  } else {
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

    boot_matrix <- matrix(NA_real_, nrow = B, ncol = length(observed))
    colnames(boot_matrix) <- stat_names
    completed <- 0L
  }

  group_indices <- split(seq_len(nrow(X)), g)

  save_checkpoint <- function(completed_resamples) {
    if (is.null(checkpoint_file)) {
      return(invisible(NULL))
    }

    rng_state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }

    saveRDS(
      list(
        observed = observed,
        stat_names = stat_names,
        boot_matrix = boot_matrix,
        completed = completed_resamples,
        B = B,
        conf_level = conf_level,
        min_finite = min_finite,
        null_value = null_value,
        rng_state = rng_state
      ),
      checkpoint_file
    )
  }

  if (!checkpoint_exists) {
    save_checkpoint(0L)
  }

  pb <- NULL
  if (isTRUE(progress)) {
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    setTxtProgressBar(pb, completed)
  }

  start_b <- completed + 1L
  if (start_b <= B) {
    for (b in seq.int(start_b, B)) {
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

      if (!is.null(checkpoint_file) && (b %% checkpoint_every == 0L || b == B)) {
        save_checkpoint(b)
      }
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
    conf_level = conf_level,
    checkpoint_file = checkpoint_file,
    checkpoint_every = checkpoint_every,
    resumed_from_checkpoint = resumed_from_checkpoint
  )
}

build_pairwise_fdr_outputs <- function(branch_object,
                                       X,
                                       g,
                                       branch_label,
                                       selection_rule,
                                       B,
                                       bootstrap_seed,
                                       checkpoint_file,
                                       checkpoint_every,
                                       min_finite = NULL) {
  X <- as.matrix(X)
  varnames <- colnames(X)
  pair_labels <- branch_object$pairwise$pair_label
  stat_names <- c(paste0("delta_A::", pair_labels), paste0("r_R_het::", pair_labels))

  if (is.null(min_finite)) {
    min_finite <- max(50L, floor(B * 0.50))
  }

  pairwise_stat_fn <- function(Xb, gb) {
    result_b <- iwaba_full(Xb, gb)
    delta_A_b <- iwaba_pairwise_dominance_contrast(result_b, varnames = varnames)
    names(delta_A_b) <- iwaba_pair_labels(result_b$level2$pairs, varnames)

    r_R_het_b <- result_b$level5$r_R_het
    names(r_R_het_b) <- iwaba_pair_labels(result_b$level2$pairs, varnames)

    stats::setNames(
      c(delta_A_b[pair_labels], r_R_het_b[pair_labels]),
      stat_names
    )
  }

  boot <- checkpointed_within_group_bootstrap_vector(
    X = X,
    g = g,
    statistic_fn = pairwise_stat_fn,
    B = B,
    conf_level = 0.95,
    seed = bootstrap_seed,
    min_finite = min_finite,
    null_value = 0,
    checkpoint_file = checkpoint_file,
    checkpoint_every = checkpoint_every,
    progress = TRUE
  )

  pairwise_long <- boot$ci
  name_parts <- strsplit(pairwise_long$stat_name, "::", fixed = TRUE)
  pairwise_long$target <- vapply(name_parts, `[`, character(1), 1)
  pairwise_long$pair_label <- vapply(name_parts, `[`, character(1), 2)
  pairwise_long$branch <- branch_label
  pairwise_long$selection_rule <- selection_rule
  pairwise_long$bootstrap_B <- B
  pairwise_long$bootstrap_seed <- bootstrap_seed
  pairwise_long$checkpoint_file <- if (is.null(checkpoint_file)) "" else checkpoint_file
  pairwise_long$checkpoint_every <- checkpoint_every

  pair_map <- branch_object$pairwise[, c("pair_label", "var1", "var2"), drop = FALSE]
  pairwise_long <- merge(pair_map, pairwise_long, by = "pair_label", all.y = TRUE, sort = FALSE)
  pairwise_long <- pairwise_long[, c(
    "branch", "selection_rule", "target", "pair_label", "var1", "var2",
    "estimate", "ci_lower", "ci_upper", "ci_defined", "n_finite",
    "boot_defined_rate", "reject_null", "bootstrap_B", "bootstrap_seed",
    "checkpoint_file", "checkpoint_every"
  )]

  pairwise_long$p_value <- vapply(seq_len(nrow(pairwise_long)), function(i) {
    bootstrap_two_sided_p(
      samples = boot$boot_matrix[, i],
      null_value = 0,
      min_finite = min_finite
    )
  }, numeric(1))

  pairwise_long$q_value <- NA_real_
  pairwise_long$q_le_0_05 <- FALSE
  pairwise_long$q_le_0_10 <- FALSE

  for (target_name in unique(pairwise_long$target)) {
    idx <- pairwise_long$target == target_name
    pairwise_long$q_value[idx] <- p.adjust(pairwise_long$p_value[idx], method = "BH")
    pairwise_long$q_le_0_05[idx] <- !is.na(pairwise_long$q_value[idx]) & pairwise_long$q_value[idx] <= 0.05
    pairwise_long$q_le_0_10[idx] <- !is.na(pairwise_long$q_value[idx]) & pairwise_long$q_value[idx] <= 0.10
  }

  observed_pairwise <- branch_object$pairwise[, c(
    "pair_label", "var1", "var2", "delta_A",
    "dominance_waba", "dominance_iwaba", "flipped"
  ), drop = FALSE]
  observed_level5 <- branch_object$level5[, c(
    "pair_label", "r_R_billard", "r_R_het", "r_R_het_defined"
  ), drop = FALSE]
  observed_context <- merge(observed_pairwise, observed_level5, by = "pair_label", all = TRUE, sort = FALSE)

  delta_df <- pairwise_long[pairwise_long$target == "delta_A", c(
    "pair_label", "estimate", "ci_lower", "ci_upper", "ci_defined",
    "n_finite", "boot_defined_rate", "reject_null", "p_value", "q_value",
    "q_le_0_05", "q_le_0_10"
  ), drop = FALSE]
  names(delta_df)[names(delta_df) != "pair_label"] <- paste0("delta_A_", names(delta_df)[names(delta_df) != "pair_label"])

  rhet_df <- pairwise_long[pairwise_long$target == "r_R_het", c(
    "pair_label", "estimate", "ci_lower", "ci_upper", "ci_defined",
    "n_finite", "boot_defined_rate", "reject_null", "p_value", "q_value",
    "q_le_0_05", "q_le_0_10"
  ), drop = FALSE]
  names(rhet_df)[names(rhet_df) != "pair_label"] <- paste0("r_R_het_", names(rhet_df)[names(rhet_df) != "pair_label"])

  pairwise_wide <- merge(observed_context, delta_df, by = "pair_label", all.x = TRUE, sort = FALSE)
  pairwise_wide <- merge(pairwise_wide, rhet_df, by = "pair_label", all.x = TRUE, sort = FALSE)
  pairwise_wide$branch <- branch_label
  pairwise_wide$selection_rule <- selection_rule
  pairwise_wide$bootstrap_B <- B
  pairwise_wide$bootstrap_seed <- bootstrap_seed
  pairwise_wide$checkpoint_file <- if (is.null(checkpoint_file)) "" else checkpoint_file
  pairwise_wide$checkpoint_every <- checkpoint_every
  pairwise_wide <- pairwise_wide[, c(
    "branch", "selection_rule", "pair_label", "var1", "var2",
    "delta_A", "dominance_waba", "dominance_iwaba", "flipped",
    "r_R_billard", "r_R_het", "r_R_het_defined",
    "delta_A_estimate", "delta_A_ci_lower", "delta_A_ci_upper", "delta_A_ci_defined",
    "delta_A_n_finite", "delta_A_boot_defined_rate", "delta_A_reject_null",
    "delta_A_p_value", "delta_A_q_value", "delta_A_q_le_0_05", "delta_A_q_le_0_10",
    "r_R_het_estimate", "r_R_het_ci_lower", "r_R_het_ci_upper", "r_R_het_ci_defined",
    "r_R_het_n_finite", "r_R_het_boot_defined_rate", "r_R_het_reject_null",
    "r_R_het_p_value", "r_R_het_q_value", "r_R_het_q_le_0_05", "r_R_het_q_le_0_10",
    "bootstrap_B", "bootstrap_seed", "checkpoint_file", "checkpoint_every"
  )]

  pairwise_summary <- data.frame(
    branch = branch_label,
    selection_rule = selection_rule,
    n_pairs = nrow(branch_object$pairwise),
    bootstrap_B = B,
    bootstrap_seed = bootstrap_seed,
    checkpoint_file = if (is.null(checkpoint_file)) "" else checkpoint_file,
    checkpoint_every = checkpoint_every,
    resumed_from_checkpoint = boot$resumed_from_checkpoint,
    delta_A_q_le_0_05 = sum(pairwise_wide$delta_A_q_le_0_05, na.rm = TRUE),
    delta_A_q_le_0_10 = sum(pairwise_wide$delta_A_q_le_0_10, na.rm = TRUE),
    r_R_het_q_le_0_05 = sum(pairwise_wide$r_R_het_q_le_0_05, na.rm = TRUE),
    r_R_het_q_le_0_10 = sum(pairwise_wide$r_R_het_q_le_0_10, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  list(
    detail_wide = pairwise_wide,
    detail_long = pairwise_long,
    summary = pairwise_summary,
    checkpoint_file = checkpoint_file,
    resumed_from_checkpoint = boot$resumed_from_checkpoint
  )
}

write_top50_figure <- function(branch, figure_filename, output_dir, tag = NULL) {
  l1 <- branch$result$level1
  l4 <- branch$result$level4
  pdf_path <- iwaba_open_pdf_device(figure_filename, output_dir = output_dir, tag = tag, width = 14, height = 12)

  par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))

  plot(l1$eta_results$eta_B_waba, l1$eta_results$eta_B_iwaba,
       xlab = expression(eta[B]^"(WABA)"),
       ylab = expression(eta[B]^"(I-WABA)"),
       main = "(a) Level I: Top-50 Entity Analysis",
       pch = 19, col = "steelblue", cex = 1.6,
       xlim = c(0, 1), ylim = c(0, 1))
  abline(0, 1, col = "red", lty = 2, lwd = 2.5)
  grid(col = "gray80", lty = 1)

  hist(l1$eta_results$delta_eta, breaks = 15, col = "#DD8452",
       border = "darkgray",
       xlim = c(min(0, min(l1$eta_results$delta_eta)), max(l1$eta_results$delta_eta) * 1.05),
       main = expression("(b) Distribution of " * Delta * eta[B] * " (Top-50)"),
       xlab = expression(Delta * eta[B]),
       ylab = "Frequency")
  abline(v = mean(l1$eta_results$delta_eta), col = "blue", lwd = 2.5)
  abline(v = 0, col = "red", lty = 2, lwd = 2)
  legend("topright", legend = c("Mean", "Zero"), col = c("blue", "red"), lty = c(1, 2), bty = "n")

  mean_radius_top50 <- rowMeans(branch$result$gq$radii)
  sorted_idx <- order(mean_radius_top50, decreasing = TRUE)
  sorted_names <- names(sort(mean_radius_top50, decreasing = TRUE))
  boxplot(t(branch$result$gq$radii[sorted_idx, , drop = FALSE]),
          las = 2, col = "#55A868",
          main = "(c) Within-Group SD by Cancer Type",
          ylab = "Within-group SD",
          cex.axis = 0.75,
          names = sorted_names,
          border = "darkgray", pch = 16)
  grid(col = "gray80", lty = 1)

  pi_R_sorted <- sort(l4$pi_R, decreasing = TRUE)
  barplot(pi_R_sorted,
          col = "#4C72B0",
          main = expression("(d) Level IV: Radius Share " * pi[R] * " (Top-50)"),
          ylab = expression(pi[R]),
          xlab = "Top-50 probes (sorted)",
          border = "darkgray",
          ylim = c(0, max(pi_R_sorted) * 1.1))
  abline(h = mean(l4$pi_R), col = "red", lty = 2, lwd = 2)
  legend("topright", legend = "Mean", col = "red", lty = 2, bty = "n")
  grid(col = "gray80", lty = 1)

  dev.off()
  pdf_path
}

write_selection_sensitivity_figure <- function(comparison_df,
                                               figure_filename,
                                               output_dir,
                                               tag = NULL) {
  pdf_path <- iwaba_open_pdf_device(figure_filename, output_dir = output_dir, tag = tag, width = 14, height = 10)

  par(mfrow = c(2, 2), mar = c(8, 5, 4, 2))
  branches <- comparison_df$branch

  barplot(comparison_df$mean_eta_B_int,
          names.arg = branches,
          las = 2,
          col = "#4C72B0",
          border = "darkgray",
          main = expression("(a) Mean " * eta[B]^"(Int)"),
          ylab = expression(Mean ~ eta[B]^"(Int)"))
  grid(col = "gray80", lty = 1)

  barplot(comparison_df$mean_E_R_het_int,
          names.arg = branches,
          las = 2,
          col = "#55A868",
          border = "darkgray",
          main = expression("(b) Mean " * E[R * "," * het]^"(Int)"),
          ylab = expression(Mean ~ E[R * "," * het]^"(Int)"))
  grid(col = "gray80", lty = 1)

  selected_mask <- !is.na(comparison_df$mean_r_R_het_signed)
  barplot(comparison_df$mean_r_R_het_signed[selected_mask],
          names.arg = branches[selected_mask],
          las = 2,
          col = "#C44E52",
          border = "darkgray",
          main = expression("(c) Mean signed " * r[R * "," * het]),
          ylab = expression(Mean ~ signed ~ r[R * "," * het]))
  abline(h = 0, col = "gray40", lwd = 1)
  grid(col = "gray80", lty = 1)

  barplot(comparison_df$bf_bh_prop,
          names.arg = branches,
          las = 2,
          col = "#8172B2",
          border = "darkgray",
          main = "(d) BF-BH discovery proportion",
          ylab = "Proportion BH-significant")
  grid(col = "gray80", lty = 1)

  dev.off()
  pdf_path
}

args <- commandArgs(trailingOnly = TRUE)
BOOT_B <- as.integer(parse_arg_value(args, "--B=", 399L))
ALPHA <- as.numeric(parse_arg_value(args, "--alpha=", 0.05))
BOOT_SEED <- as.integer(parse_arg_value(args, "--seed=", 20260321L))
MAX_PAIR_BOOT <- as.integer(parse_arg_value(args, "--max-pairs=", 10L))
CHECKPOINT_EVERY <- as.integer(parse_arg_value(args, "--checkpoint-every=", 25L))
OUTPUT_DIR <- parse_arg_value(args, "--output-dir=", "results")
RUN_TAG <- iwaba_normalize_tag(parse_arg_value(args, "--tag=", ""))
RUN_STAGE <- tolower(parse_arg_value(args, "--stage=", "full"))

VALID_STAGES <- c("full", "top50", "top25", "moderate50", "allgenes", "pairwise_fdr")
if (!RUN_STAGE %in% VALID_STAGES) {
  stop("`--stage=` must be one of: full, top50, top25, moderate50, allgenes, pairwise_fdr.")
}

RUN_TOP50 <- RUN_STAGE %in% c("full", "top50")
RUN_TOP25 <- RUN_STAGE %in% c("full", "top25")
RUN_MODERATE50 <- RUN_STAGE %in% c("full", "moderate50")
RUN_ALLGENES <- RUN_STAGE %in% c("full", "allgenes")
RUN_COMPARISON <- RUN_STAGE == "full"
RUN_PAIRCOUNT <- RUN_STAGE %in% c("full", "allgenes")
RUN_PAIRWISE_FDR <- RUN_STAGE == "pairwise_fdr"

TOP50_BOOT_SEED <- BOOT_SEED
TOP25_BOOT_SEED <- BOOT_SEED + 1L
MODERATE50_SELECTION_SEED <- BOOT_SEED + 2L
MODERATE50_BOOT_SEED <- BOOT_SEED + 3L
PAIRWISE_TOP50_BOOT_SEED <- BOOT_SEED
PAIRWISE_TOP25_BOOT_SEED <- BOOT_SEED + 1L
PAIRWISE_ILLUSTRATION_BOOT_SEED <- BOOT_SEED + 2L

TOP50_FILE <- "cancer14_top50_v3.csv"
ALLGENES_FILE <- "cancer14_allgenes_v3.csv"
TOP25_FILE <- "cancer14_top25_v1.csv"
MODERATE50_FILE <- "cancer14_moderate50_v1.csv"
MODERATE50_SELECTED_FILE <- "cancer14_selected_genes_moderate50_v1.csv"
SELECTION_COMPARISON_FILE <- "cancer14_selection_sensitivity_v1.csv"
ALLGENES_PAIRCOUNT_FILE <- "cancer14_allgenes_paircount_v1.csv"
PAIRWISE_TOP50_FILE <- "cancer14_pairwise_fdr_top50_v1.csv"
PAIRWISE_TOP50_SUMMARY_FILE <- "cancer14_pairwise_fdr_top50_summary_v1.csv"
PAIRWISE_TOP25_FILE <- "cancer14_pairwise_fdr_top25_v1.csv"
ALLGENES_PAIRWISE_ILLUSTRATION_FILE <- "cancer14_allgenes_pairwise_illustration_v1.csv"
PAIRWISE_TOP50_CHECKPOINT_FILE <- "cancer14_pairwise_fdr_top50_checkpoint_v1.rds"
PAIRWISE_TOP25_CHECKPOINT_FILE <- "cancer14_pairwise_fdr_top25_checkpoint_v1.rds"
ALLGENES_PAIRWISE_ILLUSTRATION_CHECKPOINT_FILE <- "cancer14_allgenes_pairwise_illustration_checkpoint_v1.rds"
METADATA_FILE <- "cancer14_metadata_v1.csv"
TOP50_FIGURE_FILE <- "fig_14cancer_top50_v3.pdf"
SELECTION_FIGURE_FILE <- "fig_14cancer_selection_sensitivity_v1.pdf"
RUNLOG_FILE <- file.path(OUTPUT_DIR, "realdata_runlog_v2.md")

SCRIPT_INVOCATION <- paste(c("real_data_14cancer.R", args), collapse = " ")
if (!nzchar(SCRIPT_INVOCATION)) {
  SCRIPT_INVOCATION <- "real_data_14cancer.R"
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
  top50_summary = NULL,
  top25_summary = NULL,
  moderate50_summary = NULL,
  allgenes_summary = NULL,
  allgenes_paircount = NULL,
  pairwise_top50_summary = NULL,
  pairwise_top25_summary = NULL,
  pairwise_illustration_summary = NULL
)

cat("================================================================\n")
cat("Loading 14-Cancer data (Ramaswamy et al., 2001, PNAS)\n")
cat("================================================================\n")

xtrain <- as.matrix(read.table("data/14cancer.xtrain"))
xtest <- as.matrix(read.table("data/14cancer.xtest"))
ytrain <- scan("data/14cancer.ytrain", what = integer(), quiet = TRUE)
ytest <- scan("data/14cancer.ytest", what = integer(), quiet = TRUE)

X_all <- t(cbind(xtrain, xtest))
y_all <- c(ytrain, ytest)

CANCER_NAMES <- c(
  "Breast", "Prostate", "Lung", "Colorectal", "Lymphoma",
  "Bladder", "Melanoma", "Uterus", "Leukemia", "Renal",
  "Pancreas", "Ovary", "Mesothelioma", "CNS"
)

cancer_labels <- factor(CANCER_NAMES[y_all], levels = CANCER_NAMES)
probe_ids <- paste0("Probe_", seq_len(ncol(X_all)))
colnames(X_all) <- probe_ids

cat(sprintf("  Raw combined data: n = %d, p = %d, K = %d\n",
            nrow(X_all), ncol(X_all), nlevels(cancer_labels)))
cat(sprintf("  Stage = %s, B = %d, alpha = %.3f, seed = %d, max pair boot = %d, checkpoint every = %d\n",
            RUN_STAGE, BOOT_B, ALPHA, BOOT_SEED, MAX_PAIR_BOOT, CHECKPOINT_EVERY))

gene_var <- apply(X_all, 2, var)
X_all <- X_all[, gene_var > 1e-10, drop = FALSE]
X_scaled <- scale(X_all)
X_scaled <- as.matrix(X_scaled)
colnames(X_scaled) <- colnames(X_all)

cat(sprintf("  After variance filter and global z-score: p = %d\n", ncol(X_scaled)))

grand_mean <- colMeans(X_scaled)
BSS <- numeric(ncol(X_scaled))
WSS <- numeric(ncol(X_scaled))

for (ct in levels(cancer_labels)) {
  mask <- cancer_labels == ct
  n_k <- sum(mask)
  group_mean <- colMeans(X_scaled[mask, , drop = FALSE])
  BSS <- BSS + n_k * (group_mean - grand_mean)^2
  WSS <- WSS + colSums(sweep(X_scaled[mask, , drop = FALSE], 2, group_mean)^2)
}

F_scores <- (BSS / (nlevels(cancer_labels) - 1L)) / (WSS / (nrow(X_scaled) - nlevels(cancer_labels)))
F_scores[!is.finite(F_scores)] <- 0
names(F_scores) <- colnames(X_scaled)

rank_order <- order(F_scores, decreasing = TRUE)
rank_df <- data.frame(
  variable = colnames(X_scaled)[rank_order],
  F_score = F_scores[rank_order],
  rank = seq_along(rank_order),
  stringsAsFactors = FALSE
)

top50_rank_positions <- seq_len(50L)
top25_rank_positions <- seq_len(25L)
band_start_rank <- floor(0.30 * nrow(rank_df)) + 1L
band_end_rank <- ceiling(0.70 * nrow(rank_df))
moderate_band_positions <- seq.int(band_start_rank, band_end_rank)

set.seed(MODERATE50_SELECTION_SEED)
moderate50_rank_positions <- sort(sample(moderate_band_positions, size = 50L, replace = FALSE))

top50_genes_df <- build_selected_genes_df(rank_df, top50_rank_positions, "top50")
top25_genes_df <- build_selected_genes_df(rank_df, top25_rank_positions, "top25")
moderate50_genes_df <- build_selected_genes_df(rank_df, moderate50_rank_positions, "moderate50")
moderate50_genes_df$band_start_rank <- band_start_rank
moderate50_genes_df$band_end_rank <- band_end_rank
moderate50_genes_df$selection_seed <- MODERATE50_SELECTION_SEED
moderate50_genes_df <- moderate50_genes_df[, c(
  "branch", "variable", "F_score", "rank",
  "band_start_rank", "band_end_rank", "selection_seed"
)]

TOP50_SELECTION_RULE <- "Top 50 probes ranked by one-way ANOVA F statistic after variance filtering and global z-score normalization."
TOP25_SELECTION_RULE <- "Top 25 probes ranked by the same one-way ANOVA F statistic used for the top-50 branch."
MODERATE50_SELECTION_RULE <- sprintf(
  paste(
    "Sample 50 probes uniformly without replacement from the middle 40%% of the ANOVA F ranking",
    "(ranks %d-%d of %d) after excluding the top and bottom 30%%."
  ),
  band_start_rank, band_end_rank, nrow(rank_df)
)
ALLGENES_SELECTION_RULE <- "Use all retained probes after variance filtering and global z-score normalization; no adaptive subset selection beyond preprocessing."

update_metadata <- function() {
  metadata <- list(
    script = "real_data_14cancer.R",
    dataset = "14cancer",
    run_stage = RUN_STAGE,
    output_dir = normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE),
    run_tag = if (is.null(RUN_TAG)) "" else RUN_TAG,
    alpha = ALPHA,
    bootstrap_resamples_selected_branches = BOOT_B,
    bootstrap_max_pairwise_subset = MAX_PAIR_BOOT,
    base_seed = BOOT_SEED,
    top50_boot_seed = TOP50_BOOT_SEED,
    top25_boot_seed = TOP25_BOOT_SEED,
    moderate50_selection_seed = MODERATE50_SELECTION_SEED,
    moderate50_boot_seed = MODERATE50_BOOT_SEED,
    selection_rule_top50 = TOP50_SELECTION_RULE,
    selection_rule_top25 = TOP25_SELECTION_RULE,
    selection_rule_moderate50 = MODERATE50_SELECTION_RULE,
    selection_rule_allgenes = ALLGENES_SELECTION_RULE,
    moderate50_band_start_rank = band_start_rank,
    moderate50_band_end_rank = band_end_rank,
    moderate50_band_fraction = "middle 40%",
    multiplicity_procedure = sprintf("Brown-Forsythe variablewise screening with BH adjustment at alpha = %.3f", ALPHA),
    pairwise_scope_top50 = if (is.null(analysis_state$pairwise_top50_summary)) {
      sprintf("Restricted within-group bootstrap over up to %d selected pairs only", MAX_PAIR_BOOT)
    } else {
      sprintf(
        "Exhaustive within-group bootstrap over all 1,225 top-50 pairs with BH adjustment; checkpoints every %d resamples",
        analysis_state$pairwise_top50_summary$checkpoint_every
      )
    },
    pairwise_scope_top25 = if (is.null(analysis_state$pairwise_top25_summary)) {
      sprintf("Restricted within-group bootstrap over up to %d selected pairs only", MAX_PAIR_BOOT)
    } else {
      sprintf(
        "Exhaustive within-group bootstrap over all 300 top-25 pairs with BH adjustment; checkpoints every %d resamples",
        analysis_state$pairwise_top25_summary$checkpoint_every
      )
    },
    pairwise_scope_moderate50 = sprintf("Restricted within-group bootstrap over up to %d selected pairs only", MAX_PAIR_BOOT),
    pairwise_scope_allgenes = if (is.null(analysis_state$pairwise_illustration_summary)) {
      "No exhaustive all-gene pairwise bootstrap inference in this milestone; only BF-BH heterogeneous-gene count and choose(m, 2) candidate-pair scope are saved."
    } else {
      "No exhaustive all-gene pairwise bootstrap inference; a separate illustration uses the top 50 BH-significant genes ranked by E_R_het^(Int)."
    },
    started_at = analysis_state$started_at,
    metadata_updated_at = timestamp_now(),
    output_timestamps = format_output_registry(output_registry),
    outputs_created = if (nrow(output_registry) == 0L) "" else paste(output_registry$path, collapse = "; ")
  )

  if (!is.null(analysis_state$top50_summary)) {
    metadata$top50_G <- analysis_state$top50_summary$G
    metadata$top50_radius_prop <- analysis_state$top50_summary$radius_prop
    metadata$top50_mean_eta_B_int <- analysis_state$top50_summary$mean_eta_B_int
    metadata$top50_mean_E_R_het_int <- analysis_state$top50_summary$mean_E_R_het_int
    metadata$top50_bf_bh_count <- analysis_state$top50_summary$bf_bh_count
  }

  if (!is.null(analysis_state$top25_summary)) {
    metadata$top25_G <- analysis_state$top25_summary$G
    metadata$top25_radius_prop <- analysis_state$top25_summary$radius_prop
    metadata$top25_mean_eta_B_int <- analysis_state$top25_summary$mean_eta_B_int
    metadata$top25_mean_E_R_het_int <- analysis_state$top25_summary$mean_E_R_het_int
    metadata$top25_bf_bh_count <- analysis_state$top25_summary$bf_bh_count
  }

  if (!is.null(analysis_state$moderate50_summary)) {
    metadata$moderate50_G <- analysis_state$moderate50_summary$G
    metadata$moderate50_radius_prop <- analysis_state$moderate50_summary$radius_prop
    metadata$moderate50_mean_eta_B_int <- analysis_state$moderate50_summary$mean_eta_B_int
    metadata$moderate50_mean_E_R_het_int <- analysis_state$moderate50_summary$mean_E_R_het_int
    metadata$moderate50_bf_bh_count <- analysis_state$moderate50_summary$bf_bh_count
  }

  if (!is.null(analysis_state$allgenes_summary)) {
    metadata$allgenes_G <- analysis_state$allgenes_summary$G
    metadata$allgenes_radius_prop <- analysis_state$allgenes_summary$radius_prop
    metadata$allgenes_mean_eta_B_int <- analysis_state$allgenes_summary$mean_eta_B_int
    metadata$allgenes_mean_E_R_het_int <- analysis_state$allgenes_summary$mean_E_R_het_int
    metadata$allgenes_bf_bh_count <- analysis_state$allgenes_summary$bf_bh_count
  }

  if (!is.null(analysis_state$allgenes_paircount)) {
    metadata$allgenes_bh_positive_gene_count <- analysis_state$allgenes_paircount$bh_positive_gene_count
    metadata$allgenes_candidate_pair_count <- analysis_state$allgenes_paircount$candidate_pair_count
  }

  if (!is.null(analysis_state$pairwise_top50_summary)) {
    metadata$pairwise_fdr_top50_bootstrap_B <- analysis_state$pairwise_top50_summary$bootstrap_B
    metadata$pairwise_fdr_top50_bootstrap_seed <- analysis_state$pairwise_top50_summary$bootstrap_seed
    metadata$pairwise_fdr_top50_checkpoint_every <- analysis_state$pairwise_top50_summary$checkpoint_every
    metadata$pairwise_fdr_top50_delta_A_q_le_0_05 <- analysis_state$pairwise_top50_summary$delta_A_q_le_0_05
    metadata$pairwise_fdr_top50_delta_A_q_le_0_10 <- analysis_state$pairwise_top50_summary$delta_A_q_le_0_10
    metadata$pairwise_fdr_top50_r_R_het_q_le_0_05 <- analysis_state$pairwise_top50_summary$r_R_het_q_le_0_05
    metadata$pairwise_fdr_top50_r_R_het_q_le_0_10 <- analysis_state$pairwise_top50_summary$r_R_het_q_le_0_10
  }

  if (!is.null(analysis_state$pairwise_top25_summary)) {
    metadata$pairwise_fdr_top25_bootstrap_B <- analysis_state$pairwise_top25_summary$bootstrap_B
    metadata$pairwise_fdr_top25_bootstrap_seed <- analysis_state$pairwise_top25_summary$bootstrap_seed
    metadata$pairwise_fdr_top25_delta_A_q_le_0_05 <- analysis_state$pairwise_top25_summary$delta_A_q_le_0_05
    metadata$pairwise_fdr_top25_delta_A_q_le_0_10 <- analysis_state$pairwise_top25_summary$delta_A_q_le_0_10
    metadata$pairwise_fdr_top25_r_R_het_q_le_0_05 <- analysis_state$pairwise_top25_summary$r_R_het_q_le_0_05
    metadata$pairwise_fdr_top25_r_R_het_q_le_0_10 <- analysis_state$pairwise_top25_summary$r_R_het_q_le_0_10
  }

  if (!is.null(analysis_state$pairwise_illustration_summary)) {
    metadata$pairwise_illustration_bootstrap_B <- analysis_state$pairwise_illustration_summary$bootstrap_B
    metadata$pairwise_illustration_bootstrap_seed <- analysis_state$pairwise_illustration_summary$bootstrap_seed
    metadata$pairwise_illustration_delta_A_q_le_0_05 <- analysis_state$pairwise_illustration_summary$delta_A_q_le_0_05
    metadata$pairwise_illustration_delta_A_q_le_0_10 <- analysis_state$pairwise_illustration_summary$delta_A_q_le_0_10
    metadata$pairwise_illustration_r_R_het_q_le_0_05 <- analysis_state$pairwise_illustration_summary$r_R_het_q_le_0_05
    metadata$pairwise_illustration_r_R_het_q_le_0_10 <- analysis_state$pairwise_illustration_summary$r_R_het_q_le_0_10
  }

  metadata_paths <- iwaba_save_metadata(
    metadata = metadata,
    metadata_filename = METADATA_FILE,
    output_dir = OUTPUT_DIR,
    tag = RUN_TAG
  )

  output_registry <<- output_registry[output_registry$label != "metadata", , drop = FALSE]
  output_registry <<- register_dual_output(
    output_registry,
    "metadata",
    list(csv = metadata_paths$csv, rds = metadata_paths$rds, timestamp = timestamp_now())
  )
  metadata_paths
}

top50_branch <- NULL
top25_branch <- NULL
moderate50_branch <- NULL
allgenes_branch <- NULL

if (RUN_PAIRWISE_FDR) {
  cat("\n================================================================\n")
  cat("Pairwise FDR milestone: loading saved branch objects\n")
  cat("================================================================\n")

  top50_source_path <- iwaba_output_path(sub("\\.csv$", ".rds", TOP50_FILE), output_dir = OUTPUT_DIR, tag = RUN_TAG, create_dir = FALSE)
  top25_source_path <- iwaba_output_path(sub("\\.csv$", ".rds", TOP25_FILE), output_dir = OUTPUT_DIR, tag = RUN_TAG, create_dir = FALSE)
  allgenes_source_path <- iwaba_output_path(sub("\\.csv$", ".rds", ALLGENES_FILE), output_dir = OUTPUT_DIR, tag = RUN_TAG, create_dir = FALSE)

  required_paths <- c(top50_source_path, top25_source_path, allgenes_source_path)
  if (!all(file.exists(required_paths))) {
    stop(
      paste(
        "Pairwise FDR stage requires saved branch objects from the earlier 14-Cancer milestone.",
        "Missing:",
        paste(required_paths[!file.exists(required_paths)], collapse = "; ")
      )
    )
  }

  top50_branch <- readRDS(top50_source_path)
  top25_branch <- readRDS(top25_source_path)
  allgenes_branch <- readRDS(allgenes_source_path)

  analysis_state$top50_summary <- top50_branch$summary
  analysis_state$top25_summary <- top25_branch$summary
  analysis_state$allgenes_summary <- allgenes_branch$summary

  cat(sprintf("  Loaded top-50 branch from %s\n", top50_source_path))
  cat(sprintf("  Loaded top-25 branch from %s\n", top25_source_path))
  cat(sprintf("  Loaded all-gene branch from %s\n", allgenes_source_path))

  cat("\n================================================================\n")
  cat("Task A: Exhaustive top-50 pairwise BH follow-up\n")
  cat("================================================================\n")

  X_top50_pairwise <- X_scaled[, match(top50_branch$selected_genes$variable, colnames(X_scaled)), drop = FALSE]
  colnames(X_top50_pairwise) <- top50_branch$selected_genes$variable
  top50_checkpoint_path <- iwaba_output_path(PAIRWISE_TOP50_CHECKPOINT_FILE, output_dir = OUTPUT_DIR, tag = RUN_TAG)

  top50_pairwise_fdr <- build_pairwise_fdr_outputs(
    branch_object = top50_branch,
    X = X_top50_pairwise,
    g = cancer_labels,
    branch_label = "top50",
    selection_rule = top50_branch$selection_rule,
    B = BOOT_B,
    bootstrap_seed = PAIRWISE_TOP50_BOOT_SEED,
    checkpoint_file = top50_checkpoint_path,
    checkpoint_every = CHECKPOINT_EVERY
  )
  analysis_state$pairwise_top50_summary <- top50_pairwise_fdr$summary

  top50_pairwise_paths <- save_dual_output(
    csv_object = top50_pairwise_fdr$detail_wide,
    rds_object = top50_pairwise_fdr,
    csv_filename = PAIRWISE_TOP50_FILE,
    output_dir = OUTPUT_DIR,
    tag = RUN_TAG
  )
  output_registry <- register_dual_output(output_registry, "pairwise_fdr_top50", top50_pairwise_paths)
  output_registry <- register_single_output(output_registry, "pairwise_fdr_top50_checkpoint", top50_checkpoint_path, type = "checkpoint_rds")

  top50_pairwise_summary_paths <- save_dual_output(
    csv_object = top50_pairwise_fdr$summary,
    rds_object = top50_pairwise_fdr$summary,
    csv_filename = PAIRWISE_TOP50_SUMMARY_FILE,
    output_dir = OUTPUT_DIR,
    tag = RUN_TAG
  )
  output_registry <- register_dual_output(output_registry, "pairwise_fdr_top50_summary", top50_pairwise_summary_paths)

  metadata_paths <- update_metadata()

  append_runlog_entry(
    runlog_file = RUNLOG_FILE,
    subtask = "Pairwise Task A - Top-50 exhaustive BH",
    files_changed = c(
      basename(top50_pairwise_paths$csv),
      basename(top50_pairwise_paths$rds),
      basename(top50_pairwise_summary_paths$csv),
      basename(top50_pairwise_summary_paths$rds),
      basename(top50_checkpoint_path),
      basename(metadata_paths$csv),
      basename(metadata_paths$rds)
    ),
    commands_run = c(
      SCRIPT_INVOCATION,
      sprintf(
        "Task A: exhaustive within-group bootstrap over all 1,225 top-50 pairs with B=%d and checkpoints every %d resamples",
        BOOT_B, CHECKPOINT_EVERY
      )
    ),
    outputs_created = c(
      top50_pairwise_paths$csv,
      top50_pairwise_paths$rds,
      top50_pairwise_summary_paths$csv,
      top50_pairwise_summary_paths$rds,
      top50_checkpoint_path,
      metadata_paths$csv,
      metadata_paths$rds
    ),
    seeds_used = sprintf("Top-50 pairwise bootstrap seed = %d", PAIRWISE_TOP50_BOOT_SEED),
    what_remains = c(
      "Task B: top-25 exhaustive pairwise BH",
      "Task C: restricted all-gene illustration among top 50 E_R_het-ranked BH-significant genes"
    )
  )

  cat(sprintf(
    "  Top-50 pairwise BH counts: delta_A q<=0.05 = %d, delta_A q<=0.10 = %d, r_R_het q<=0.05 = %d, r_R_het q<=0.10 = %d\n",
    top50_pairwise_fdr$summary$delta_A_q_le_0_05,
    top50_pairwise_fdr$summary$delta_A_q_le_0_10,
    top50_pairwise_fdr$summary$r_R_het_q_le_0_05,
    top50_pairwise_fdr$summary$r_R_het_q_le_0_10
  ))

  cat("\n================================================================\n")
  cat("Task B: Exhaustive top-25 pairwise BH follow-up\n")
  cat("================================================================\n")

  X_top25_pairwise <- X_scaled[, match(top25_branch$selected_genes$variable, colnames(X_scaled)), drop = FALSE]
  colnames(X_top25_pairwise) <- top25_branch$selected_genes$variable
  top25_checkpoint_path <- iwaba_output_path(PAIRWISE_TOP25_CHECKPOINT_FILE, output_dir = OUTPUT_DIR, tag = RUN_TAG)

  top25_pairwise_fdr <- build_pairwise_fdr_outputs(
    branch_object = top25_branch,
    X = X_top25_pairwise,
    g = cancer_labels,
    branch_label = "top25",
    selection_rule = top25_branch$selection_rule,
    B = BOOT_B,
    bootstrap_seed = PAIRWISE_TOP25_BOOT_SEED,
    checkpoint_file = top25_checkpoint_path,
    checkpoint_every = CHECKPOINT_EVERY
  )
  analysis_state$pairwise_top25_summary <- top25_pairwise_fdr$summary

  top25_pairwise_paths <- save_dual_output(
    csv_object = top25_pairwise_fdr$detail_wide,
    rds_object = top25_pairwise_fdr,
    csv_filename = PAIRWISE_TOP25_FILE,
    output_dir = OUTPUT_DIR,
    tag = RUN_TAG
  )
  output_registry <- register_dual_output(output_registry, "pairwise_fdr_top25", top25_pairwise_paths)
  output_registry <- register_single_output(output_registry, "pairwise_fdr_top25_checkpoint", top25_checkpoint_path, type = "checkpoint_rds")

  metadata_paths <- update_metadata()

  append_runlog_entry(
    runlog_file = RUNLOG_FILE,
    subtask = "Pairwise Task B - Top-25 exhaustive BH",
    files_changed = c(
      basename(top25_pairwise_paths$csv),
      basename(top25_pairwise_paths$rds),
      basename(top25_checkpoint_path),
      basename(metadata_paths$csv),
      basename(metadata_paths$rds)
    ),
    commands_run = c(
      SCRIPT_INVOCATION,
      sprintf(
        "Task B: exhaustive within-group bootstrap over all 300 top-25 pairs with B=%d and checkpoints every %d resamples",
        BOOT_B, CHECKPOINT_EVERY
      )
    ),
    outputs_created = c(
      top25_pairwise_paths$csv,
      top25_pairwise_paths$rds,
      top25_checkpoint_path,
      metadata_paths$csv,
      metadata_paths$rds
    ),
    seeds_used = sprintf("Top-25 pairwise bootstrap seed = %d", PAIRWISE_TOP25_BOOT_SEED),
    what_remains = "Task C: restricted all-gene illustration among top 50 E_R_het-ranked BH-significant genes"
  )

  cat(sprintf(
    "  Top-25 pairwise BH counts: delta_A q<=0.05 = %d, delta_A q<=0.10 = %d, r_R_het q<=0.05 = %d, r_R_het q<=0.10 = %d\n",
    top25_pairwise_fdr$summary$delta_A_q_le_0_05,
    top25_pairwise_fdr$summary$delta_A_q_le_0_10,
    top25_pairwise_fdr$summary$r_R_het_q_le_0_05,
    top25_pairwise_fdr$summary$r_R_het_q_le_0_10
  ))

  cat("\n================================================================\n")
  cat("Task C: Restricted all-gene pairwise illustration\n")
  cat("================================================================\n")

  allgenes_screening_rank <- allgenes_branch$screening
  allgenes_screening_rank <- allgenes_screening_rank[allgenes_screening_rank$reject, , drop = FALSE]
  allgenes_screening_rank <- allgenes_screening_rank[
    order(-allgenes_screening_rank$E_R_het_int, -allgenes_screening_rank$F_score, allgenes_screening_rank$variable),
    ,
    drop = FALSE
  ]
  illustration_selected_genes <- allgenes_screening_rank[seq_len(50L), c(
    "variable", "E_R_het_int", "F_score", "p_value", "p_adj"
  ), drop = FALSE]
  illustration_selected_genes$selection_rank <- seq_len(nrow(illustration_selected_genes))

  illustration_selection_rule <- paste(
    "Top 50 all-gene probes among BH-significant heterogeneous genes, ranked by E_R_het^(Int)",
    "from the all-gene screening branch; used only as a prespecified illustration."
  )

  X_illustration <- X_scaled[, match(illustration_selected_genes$variable, colnames(X_scaled)), drop = FALSE]
  colnames(X_illustration) <- illustration_selected_genes$variable

  illustration_branch <- run_selected_branch(
    X = X_illustration,
    g = cancer_labels,
    branch_label = "allgenes_illustration_top50_erhet",
    selection_rule = illustration_selection_rule,
    f_scores_lookup = F_scores,
    alpha = ALPHA,
    B = 0L,
    bootstrap_seed = NA_integer_,
    max_pair_boot = 0L
  )
  illustration_branch$selected_genes <- illustration_selected_genes

  illustration_checkpoint_path <- iwaba_output_path(
    ALLGENES_PAIRWISE_ILLUSTRATION_CHECKPOINT_FILE,
    output_dir = OUTPUT_DIR,
    tag = RUN_TAG
  )

  illustration_pairwise_fdr <- build_pairwise_fdr_outputs(
    branch_object = illustration_branch,
    X = X_illustration,
    g = cancer_labels,
    branch_label = "allgenes_illustration_top50_erhet",
    selection_rule = illustration_selection_rule,
    B = BOOT_B,
    bootstrap_seed = PAIRWISE_ILLUSTRATION_BOOT_SEED,
    checkpoint_file = illustration_checkpoint_path,
    checkpoint_every = CHECKPOINT_EVERY
  )
  analysis_state$pairwise_illustration_summary <- illustration_pairwise_fdr$summary

  illustration_pairwise_paths <- save_dual_output(
    csv_object = illustration_pairwise_fdr$detail_wide,
    rds_object = list(
      detail_wide = illustration_pairwise_fdr$detail_wide,
      detail_long = illustration_pairwise_fdr$detail_long,
      summary = illustration_pairwise_fdr$summary,
      selected_genes = illustration_selected_genes,
      branch = illustration_branch,
      checkpoint_file = illustration_checkpoint_path
    ),
    csv_filename = ALLGENES_PAIRWISE_ILLUSTRATION_FILE,
    output_dir = OUTPUT_DIR,
    tag = RUN_TAG
  )
  output_registry <- register_dual_output(output_registry, "allgenes_pairwise_illustration", illustration_pairwise_paths)
  output_registry <- register_single_output(output_registry, "allgenes_pairwise_illustration_checkpoint", illustration_checkpoint_path, type = "checkpoint_rds")

  metadata_paths <- update_metadata()

  append_runlog_entry(
    runlog_file = RUNLOG_FILE,
    subtask = "Pairwise Task C - All-gene illustration",
    files_changed = c(
      basename(illustration_pairwise_paths$csv),
      basename(illustration_pairwise_paths$rds),
      basename(illustration_checkpoint_path),
      basename(metadata_paths$csv),
      basename(metadata_paths$rds)
    ),
    commands_run = c(
      SCRIPT_INVOCATION,
      sprintf(
        "Task C: pairwise bootstrap over the top 50 BH-significant all-gene probes ranked by E_R_het^(Int) with B=%d and checkpoints every %d resamples",
        BOOT_B, CHECKPOINT_EVERY
      )
    ),
    outputs_created = c(
      illustration_pairwise_paths$csv,
      illustration_pairwise_paths$rds,
      illustration_checkpoint_path,
      metadata_paths$csv,
      metadata_paths$rds
    ),
    seeds_used = c(
      "All-gene illustration selection is deterministic from saved screening outputs",
      sprintf("All-gene illustration bootstrap seed = %d", PAIRWISE_ILLUSTRATION_BOOT_SEED)
    ),
    what_remains = "Pairwise multiplicity-control milestone complete"
  )

  cat(sprintf(
    "  Illustration pairwise BH counts: delta_A q<=0.05 = %d, delta_A q<=0.10 = %d, r_R_het q<=0.05 = %d, r_R_het q<=0.10 = %d\n",
    illustration_pairwise_fdr$summary$delta_A_q_le_0_05,
    illustration_pairwise_fdr$summary$delta_A_q_le_0_10,
    illustration_pairwise_fdr$summary$r_R_het_q_le_0_05,
    illustration_pairwise_fdr$summary$r_R_het_q_le_0_10
  ))
}

if (RUN_TOP50) {
  cat("\n================================================================\n")
  cat("Task A (part 1): Recomputing top-50 branch\n")
  cat("================================================================\n")

  top50_idx <- match(top50_genes_df$variable, colnames(X_scaled))
  X_top50 <- X_scaled[, top50_idx, drop = FALSE]
  top50_branch <- run_selected_branch(
    X = X_top50,
    g = cancer_labels,
    branch_label = "top50",
    selection_rule = TOP50_SELECTION_RULE,
    f_scores_lookup = F_scores,
    alpha = ALPHA,
    B = BOOT_B,
    bootstrap_seed = TOP50_BOOT_SEED,
    max_pair_boot = MAX_PAIR_BOOT
  )
  top50_branch$selected_genes <- top50_genes_df
  analysis_state$top50_summary <- top50_branch$summary

  top50_paths <- save_dual_output(
    csv_object = top50_branch$export,
    rds_object = top50_branch,
    csv_filename = TOP50_FILE,
    output_dir = OUTPUT_DIR,
    tag = RUN_TAG
  )
  output_registry <- register_dual_output(output_registry, "top50", top50_paths)

  top50_fig_path <- write_top50_figure(top50_branch, TOP50_FIGURE_FILE, OUTPUT_DIR, RUN_TAG)
  output_registry <- register_pdf_output(output_registry, "top50_figure", top50_fig_path)

  cat(sprintf("  Top-50: G = %.3f, radius share = %.3f, BF-BH discoveries = %d/%d\n",
              top50_branch$summary$G,
              top50_branch$summary$radius_prop,
              top50_branch$summary$bf_bh_count,
              top50_branch$summary$p))
}

if (RUN_ALLGENES) {
  cat("\n================================================================\n")
  cat("Task A (part 2): Recomputing all-gene screening branch\n")
  cat("================================================================\n")

  allgenes_branch <- run_allgenes_branch(
    X = X_scaled,
    g = cancer_labels,
    f_scores_lookup = F_scores,
    alpha = ALPHA,
    selection_rule = ALLGENES_SELECTION_RULE
  )
  analysis_state$allgenes_summary <- allgenes_branch$summary

  allgenes_paths <- save_dual_output(
    csv_object = allgenes_branch$export,
    rds_object = allgenes_branch,
    csv_filename = ALLGENES_FILE,
    output_dir = OUTPUT_DIR,
    tag = RUN_TAG
  )
  output_registry <- register_dual_output(output_registry, "allgenes", allgenes_paths)

  cat(sprintf("  All genes: G = %.3f, radius share = %.3f, BF-BH discoveries = %d/%d\n",
              allgenes_branch$summary$G,
              allgenes_branch$summary$radius_prop,
              allgenes_branch$summary$bf_bh_count,
              allgenes_branch$summary$p))
}

if (RUN_STAGE == "full") {
  metadata_paths <- update_metadata()
  append_runlog_entry(
    runlog_file = RUNLOG_FILE,
    subtask = "Task A - Existing top-50 and all-gene branches",
    files_changed = c(
      basename(top50_paths$csv),
      basename(top50_paths$rds),
      basename(allgenes_paths$csv),
      basename(allgenes_paths$rds),
      basename(top50_fig_path),
      basename(metadata_paths$csv),
      basename(metadata_paths$rds)
    ),
    commands_run = c(
      SCRIPT_INVOCATION,
      sprintf("Task A top-50: top 50 ANOVA F-ranked probes with restricted bootstrap B=%d", BOOT_B),
      "Task A all-genes: all retained probes with Brown-Forsythe screening only"
    ),
    outputs_created = c(
      top50_paths$csv, top50_paths$rds,
      allgenes_paths$csv, allgenes_paths$rds,
      top50_fig_path,
      metadata_paths$csv, metadata_paths$rds
    ),
    seeds_used = c(
      sprintf("Top-50 bootstrap seed = %d", TOP50_BOOT_SEED),
      "All-gene branch uses no RNG"
    ),
    what_remains = c(
      "Task B: top-25 branch",
      "Task C: moderate-F 50 branch",
      "Task D: selection-sensitivity comparison",
      "Task E: all-gene candidate-pair count"
    )
  )
}

if (RUN_TOP25) {
  cat("\n================================================================\n")
  cat("Task B: Recomputing top-25 branch\n")
  cat("================================================================\n")

  top25_idx <- match(top25_genes_df$variable, colnames(X_scaled))
  X_top25 <- X_scaled[, top25_idx, drop = FALSE]
  top25_branch <- run_selected_branch(
    X = X_top25,
    g = cancer_labels,
    branch_label = "top25",
    selection_rule = TOP25_SELECTION_RULE,
    f_scores_lookup = F_scores,
    alpha = ALPHA,
    B = BOOT_B,
    bootstrap_seed = TOP25_BOOT_SEED,
    max_pair_boot = MAX_PAIR_BOOT
  )
  top25_branch$selected_genes <- top25_genes_df
  analysis_state$top25_summary <- top25_branch$summary

  top25_paths <- save_dual_output(
    csv_object = top25_branch$export,
    rds_object = top25_branch,
    csv_filename = TOP25_FILE,
    output_dir = OUTPUT_DIR,
    tag = RUN_TAG
  )
  output_registry <- register_dual_output(output_registry, "top25", top25_paths)

  metadata_paths <- update_metadata()

  append_runlog_entry(
    runlog_file = RUNLOG_FILE,
    subtask = "Task B - Top-25 branch",
    files_changed = c(
      basename(top25_paths$csv),
      basename(top25_paths$rds),
      basename(metadata_paths$csv),
      basename(metadata_paths$rds)
    ),
    commands_run = c(
      SCRIPT_INVOCATION,
      sprintf("Task B: top 25 ANOVA F-ranked probes with restricted bootstrap B=%d", BOOT_B)
    ),
    outputs_created = c(
      top25_paths$csv, top25_paths$rds,
      metadata_paths$csv, metadata_paths$rds
    ),
    seeds_used = sprintf("Top-25 bootstrap seed = %d", TOP25_BOOT_SEED),
    what_remains = c(
      "Task C: moderate-F 50 branch",
      "Task D: selection-sensitivity comparison",
      "Task E: all-gene candidate-pair count"
    )
  )

  cat(sprintf("  Top-25: G = %.3f, radius share = %.3f, BF-BH discoveries = %d/%d\n",
              top25_branch$summary$G,
              top25_branch$summary$radius_prop,
              top25_branch$summary$bf_bh_count,
              top25_branch$summary$p))
}

if (RUN_MODERATE50) {
  cat("\n================================================================\n")
  cat("Task C: Recomputing moderate-F 50 branch\n")
  cat("================================================================\n")

  moderate50_idx <- match(moderate50_genes_df$variable, colnames(X_scaled))
  X_moderate50 <- X_scaled[, moderate50_idx, drop = FALSE]
  moderate50_branch <- run_selected_branch(
    X = X_moderate50,
    g = cancer_labels,
    branch_label = "moderate50",
    selection_rule = MODERATE50_SELECTION_RULE,
    f_scores_lookup = F_scores,
    alpha = ALPHA,
    B = BOOT_B,
    bootstrap_seed = MODERATE50_BOOT_SEED,
    max_pair_boot = MAX_PAIR_BOOT
  )
  moderate50_branch$selected_genes <- moderate50_genes_df
  analysis_state$moderate50_summary <- moderate50_branch$summary

  moderate50_paths <- save_dual_output(
    csv_object = moderate50_branch$export,
    rds_object = moderate50_branch,
    csv_filename = MODERATE50_FILE,
    output_dir = OUTPUT_DIR,
    tag = RUN_TAG
  )
  output_registry <- register_dual_output(output_registry, "moderate50", moderate50_paths)

  moderate50_selected_paths <- save_dual_output(
    csv_object = moderate50_genes_df,
    rds_object = list(
      selected_genes = moderate50_genes_df,
      selection_rule = MODERATE50_SELECTION_RULE,
      selection_seed = MODERATE50_SELECTION_SEED,
      band_start_rank = band_start_rank,
      band_end_rank = band_end_rank
    ),
    csv_filename = MODERATE50_SELECTED_FILE,
    output_dir = OUTPUT_DIR,
    tag = RUN_TAG
  )
  output_registry <- register_dual_output(output_registry, "moderate50_selected_genes", moderate50_selected_paths)

  metadata_paths <- update_metadata()

  append_runlog_entry(
    runlog_file = RUNLOG_FILE,
    subtask = "Task C - Moderate-F 50 branch",
    files_changed = c(
      basename(moderate50_paths$csv),
      basename(moderate50_paths$rds),
      basename(moderate50_selected_paths$csv),
      basename(moderate50_selected_paths$rds),
      basename(metadata_paths$csv),
      basename(metadata_paths$rds)
    ),
    commands_run = c(
      SCRIPT_INVOCATION,
      sprintf(
        "Task C: sample 50 probes uniformly from middle 40%% ANOVA F ranks (%d-%d) with selection seed %d and restricted bootstrap B=%d",
        band_start_rank, band_end_rank, MODERATE50_SELECTION_SEED, BOOT_B
      )
    ),
    outputs_created = c(
      moderate50_paths$csv, moderate50_paths$rds,
      moderate50_selected_paths$csv, moderate50_selected_paths$rds,
      metadata_paths$csv, metadata_paths$rds
    ),
    seeds_used = c(
      sprintf("Moderate-F gene-selection seed = %d", MODERATE50_SELECTION_SEED),
      sprintf("Moderate-F bootstrap seed = %d", MODERATE50_BOOT_SEED)
    ),
    what_remains = c(
      "Task D: selection-sensitivity comparison",
      "Task E: all-gene candidate-pair count"
    )
  )

  cat(sprintf("  Moderate-F 50: G = %.3f, radius share = %.3f, BF-BH discoveries = %d/%d\n",
              moderate50_branch$summary$G,
              moderate50_branch$summary$radius_prop,
              moderate50_branch$summary$bf_bh_count,
              moderate50_branch$summary$p))
}

if (RUN_COMPARISON) {
  cat("\n================================================================\n")
  cat("Task D: Building selection-sensitivity comparison\n")
  cat("================================================================\n")

  comparison_df <- bind_named_rows(list(
    top50_branch$summary,
    top25_branch$summary,
    moderate50_branch$summary,
    allgenes_branch$summary
  ))
  comparison_df$conditional_on_selection <- comparison_df$branch %in% c("top50", "top25", "moderate50")
  comparison_df$comparable_pairwise_metrics <- !is.na(comparison_df$mean_abs_r_between_int)

  selection_comparison_paths <- save_dual_output(
    csv_object = comparison_df,
    rds_object = list(
      comparison = comparison_df,
      summaries = list(
        top50 = top50_branch$summary,
        top25 = top25_branch$summary,
        moderate50 = moderate50_branch$summary,
        allgenes = allgenes_branch$summary
      )
    ),
    csv_filename = SELECTION_COMPARISON_FILE,
    output_dir = OUTPUT_DIR,
    tag = RUN_TAG
  )
  output_registry <- register_dual_output(output_registry, "selection_comparison", selection_comparison_paths)

  selection_fig_path <- write_selection_sensitivity_figure(
    comparison_df = comparison_df,
    figure_filename = SELECTION_FIGURE_FILE,
    output_dir = OUTPUT_DIR,
    tag = RUN_TAG
  )
  output_registry <- register_pdf_output(output_registry, "selection_comparison_figure", selection_fig_path)

  metadata_paths <- update_metadata()

  append_runlog_entry(
    runlog_file = RUNLOG_FILE,
    subtask = "Task D - Selection-sensitivity comparison",
    files_changed = c(
      basename(selection_comparison_paths$csv),
      basename(selection_comparison_paths$rds),
      basename(selection_fig_path),
      basename(metadata_paths$csv),
      basename(metadata_paths$rds)
    ),
    commands_run = c(
      SCRIPT_INVOCATION,
      "Task D: compare top-50, top-25, moderate-F 50, and all-gene screening summaries"
    ),
    outputs_created = c(
      selection_comparison_paths$csv,
      selection_comparison_paths$rds,
      selection_fig_path,
      metadata_paths$csv,
      metadata_paths$rds
    ),
    seeds_used = "No new RNG used beyond saved branch outputs",
    what_remains = "Task E: all-gene candidate-pair count"
  )
}

if (RUN_PAIRCOUNT) {
  cat("\n================================================================\n")
  cat("Task E: Saving all-gene BF-BH pair-count scope\n")
  cat("================================================================\n")

  bh_positive_gene_count <- sum(allgenes_branch$screening$reject, na.rm = TRUE)
  candidate_pair_count <- as.numeric(bh_positive_gene_count) * as.numeric(bh_positive_gene_count - 1) / 2

  allgenes_paircount_df <- data.frame(
    branch = "allgenes",
    alpha = ALPHA,
    multiplicity_procedure = "BH on Brown-Forsythe p-values",
    bh_positive_gene_count = bh_positive_gene_count,
    candidate_pair_count = candidate_pair_count,
    stringsAsFactors = FALSE
  )
  analysis_state$allgenes_paircount <- allgenes_paircount_df

  allgenes_paircount_paths <- save_dual_output(
    csv_object = allgenes_paircount_df,
    rds_object = list(
      paircount = allgenes_paircount_df,
      screening = allgenes_branch$screening
    ),
    csv_filename = ALLGENES_PAIRCOUNT_FILE,
    output_dir = OUTPUT_DIR,
    tag = RUN_TAG
  )
  output_registry <- register_dual_output(output_registry, "allgenes_paircount", allgenes_paircount_paths)

  metadata_paths <- update_metadata()

  append_runlog_entry(
    runlog_file = RUNLOG_FILE,
    subtask = "Task E - All-gene pair-count scope",
    files_changed = c(
      basename(allgenes_paircount_paths$csv),
      basename(allgenes_paircount_paths$rds),
      basename(metadata_paths$csv),
      basename(metadata_paths$rds)
    ),
    commands_run = c(
      SCRIPT_INVOCATION,
      "Task E: count BH-positive heterogeneous genes and choose(m, 2) candidate pairs for the all-gene branch"
    ),
    outputs_created = c(
      allgenes_paircount_paths$csv,
      allgenes_paircount_paths$rds,
      metadata_paths$csv,
      metadata_paths$rds
    ),
    seeds_used = "No RNG used",
    what_remains = "14-Cancer screening/selection milestone complete"
  )

  cat(sprintf("  All-gene BF-BH positive genes = %d, candidate pairs = %.0f\n",
              bh_positive_gene_count, candidate_pair_count))
}

cat("\nSaved outputs:\n")
for (path in output_registry$path) {
  cat(sprintf("  %s\n", path))
}

cat("\nDone.\n")
