################################################################################
## I-WABA Simulation Study (v3) — Matched p-Sensitivity + Level V Patch
## Title: Interval-based Within-And-Between Analysis
## Author: Han-Ming Wu (National Chengchi University)
##
## Uses Billard's endpoint covariance (Eq. 21) and implements all five
## I-WABA levels following the workflow in Section 2.4.5.
##
## Output: simulation_results_v3.rds (full results data.frame)
##         Versioned LaTeX tables and publication-quality figures
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

OUTPUT_VERSION <- "v3"
RESULTS_RDS <- sprintf("simulation_results_%s.rds", OUTPUT_VERSION)
RESULTS_CSV <- sprintf("simulation_results_%s.csv", OUTPUT_VERSION)
FIG_INFO_GAIN_FILE <- sprintf("fig_info_gain_vs_K_%s.pdf", OUTPUT_VERSION)
FIG_ETA_IMPROVEMENT_FILE <- sprintf("fig_eta_improvement_%s.pdf", OUTPUT_VERSION)
FIG_DISPERSION_FILE <- sprintf("fig_dispersion_diagnostics_%s.pdf", OUTPUT_VERSION)
FIG_HEATMAP_FILE <- sprintf("fig_heatmap_info_gain_%s.pdf", OUTPUT_VERSION)
FIG_RADIUS_FILE <- sprintf("fig_radius_association_%s.pdf", OUTPUT_VERSION)

mean_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

format_num_or_dash <- function(x, digits = 3) {
  if (is.na(x)) {
    "--"
  } else {
    sprintf(paste0("%.", digits, "f"), x)
  }
}

format_level5_note <- function(K, prop_pos, prop_neg, prop_undef) {
  if (K == 2) {
    sprintf("K=2 sign-only: P(+)=%.2f, P(-)=%.2f, undef=%.2f",
            prop_pos, prop_neg, prop_undef)
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

# ============================================================================
# 1. Data Generation Functions
# ============================================================================

#' Heteroscedastic: different means AND different correlated within-group variances
generate_heteroscedastic <- function(n, p, K, mean_shift = 1.0,
                                      spread_corr = 0.8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_per <- rep(floor(n / K), K)
  n_per[K] <- n - sum(n_per[-K])

  group_means <- matrix(0, K, p)
  for (j in 1:p) {
    group_means[, j] <- seq(-mean_shift, mean_shift, length.out = K) * (1 + 0.2 * (j - 1))
  }

  base_sd <- seq(0.5, 2.0, length.out = K)
  Sigma_noise <- spread_corr * matrix(1, p, p) + (1 - spread_corr) * diag(p)
  group_sds <- matrix(0, K, p)
  for (k in 1:K) {
    noise <- mvrnorm(1, mu = rep(0, p), Sigma = Sigma_noise)
    group_sds[k, ] <- pmax(base_sd[k] + 0.3 * noise, 0.1)
  }

  X <- matrix(0, n, p)
  Y <- integer(n)
  idx <- 1
  for (k in 1:K) {
    for (j in 1:p) {
      X[idx:(idx + n_per[k] - 1), j] <- rnorm(n_per[k], group_means[k, j], group_sds[k, j])
    }
    Y[idx:(idx + n_per[k] - 1)] <- k
    idx <- idx + n_per[k]
  }
  list(X = X, Y = as.factor(Y))
}

#' Spread-Only: identical means, different within-group variances
generate_spread_only <- function(n, p, K, spread_range = c(0.5, 3.0),
                                  spread_corr = 0.9, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_per <- rep(floor(n / K), K)
  n_per[K] <- n - sum(n_per[-K])

  base_sd <- seq(spread_range[1], spread_range[2], length.out = K)
  Sigma_noise <- spread_corr * matrix(1, p, p) + (1 - spread_corr) * diag(p)
  group_sds <- matrix(0, K, p)
  for (k in 1:K) {
    noise <- mvrnorm(1, mu = rep(0, p), Sigma = Sigma_noise)
    group_sds[k, ] <- pmax(base_sd[k] + 0.2 * noise, 0.1)
  }

  X <- matrix(0, n, p)
  Y <- integer(n)
  idx <- 1
  for (k in 1:K) {
    for (j in 1:p) {
      X[idx:(idx + n_per[k] - 1), j] <- rnorm(n_per[k], 0, group_sds[k, j])
    }
    Y[idx:(idx + n_per[k] - 1)] <- k
    idx <- idx + n_per[k]
  }
  list(X = X, Y = as.factor(Y))
}

#' Homoscedastic: different means, identical within-group variances (control)
generate_homoscedastic <- function(n, p, K, mean_shift = 1.0,
                                    common_sd = 1.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_per <- rep(floor(n / K), K)
  n_per[K] <- n - sum(n_per[-K])

  group_means <- matrix(0, K, p)
  for (j in 1:p) {
    group_means[, j] <- seq(-mean_shift, mean_shift, length.out = K)
  }

  X <- matrix(0, n, p)
  Y <- integer(n)
  idx <- 1
  for (k in 1:K) {
    X[idx:(idx + n_per[k] - 1), ] <- mvrnorm(n_per[k], group_means[k, ],
                                                common_sd^2 * diag(p))
    Y[idx:(idx + n_per[k] - 1)] <- k
    idx <- idx + n_per[k]
  }
  list(X = X, Y = as.factor(Y))
}

#' Nonlinear Variance: spread proportional to mean (CV = 0.5)
generate_nonlinear_variance <- function(n, p, K, mean_shift = 2.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_per <- rep(floor(n / K), K)
  n_per[K] <- n - sum(n_per[-K])

  group_means <- matrix(0, K, p)
  group_sds   <- matrix(0, K, p)
  for (j in 1:p) {
    group_means[, j] <- seq(1, 1 + mean_shift, length.out = K)
    group_sds[, j]   <- 0.5 * group_means[, j]
  }

  X <- matrix(0, n, p)
  Y <- integer(n)
  idx <- 1
  for (k in 1:K) {
    for (j in 1:p) {
      X[idx:(idx + n_per[k] - 1), j] <- rnorm(n_per[k], group_means[k, j], group_sds[k, j])
    }
    Y[idx:(idx + n_per[k] - 1)] <- k
    idx <- idx + n_per[k]
  }
  list(X = X, Y = as.factor(Y))
}


# ============================================================================
# 2. Single-Replication Analysis (returns 1-row data.frame)
# ============================================================================

#' Run full I-WABA Level I–V analysis and return summary metrics
#' @param X data matrix (n x p)
#' @param Y group labels
#' @return 1-row data.frame with all diagnostic quantities
run_one_analysis <- function(X, Y) {
  result <- iwaba_full(X, Y)
  l1 <- result$level1
  l2 <- result$level2
  l3 <- result$level3
  l4 <- result$level4
  l5 <- result$level5
  eta <- l1$eta_results

  data.frame(
    # Level I
    info_gain        = l1$info_gain,
    radius_prop      = l1$radius_prop,
    mean_eta_B_waba  = mean(eta$eta_B_waba),
    mean_eta_B_int   = mean(eta$eta_B_iwaba),
    mean_delta_eta   = mean(eta$delta_eta),
    mean_E_waba      = mean(eta$E_waba),
    mean_E_int       = mean(eta$E_iwaba),
    n_reclassified   = sum(eta$class_waba != eta$class_iwaba),
    n_between_waba   = sum(eta$class_waba == "Between"),
    n_between_int    = sum(eta$class_iwaba == "Between"),
    n_equivocal_waba = sum(eta$class_waba == "Equivocal"),
    n_equivocal_int  = sum(eta$class_iwaba == "Equivocal"),
    n_within_waba    = sum(eta$class_waba == "Within"),
    n_within_int     = sum(eta$class_iwaba == "Within"),
    mean_V_C         = mean(l1$V_C),
    mean_V_R         = mean(l1$V_R),
    mean_V_W         = mean(l1$V_W),

    # Level II
    mean_abs_bw_int  = l2$mean_abs_bw_int,
    mean_abs_ww_int  = l2$mean_abs_ww_int,
    mean_abs_bw_waba = l2$mean_abs_bw_waba,
    mean_abs_ww_waba = l2$mean_abs_ww_waba,

    # Level III
    consistency_rate_int  = l3$consistency_rate_int,
    consistency_rate_waba = l3$consistency_rate_waba,

    # Level IV
    mean_pi_R        = mean(l4$pi_R),
    mean_E_C_int     = mean(l4$E_C_int),
    mean_E_R_int     = mean(l4$E_R_int),
    mean_E_R_het_int = mean(l4$E_R_het_int),
    mean_V_R_het     = mean(l4$V_R_het),

    # Level V
    mean_r_R_billard       = l5$mean_r_R_billard,
    mean_r_R_het_signed    = l5$mean_r_R_het_signed,
    prop_r_R_het_positive  = l5$prop_r_R_het_positive,
    prop_r_R_het_negative  = l5$prop_r_R_het_negative,
    prop_r_R_het_undefined = l5$prop_r_R_het_undefined,
    n_r_R_het_defined      = l5$n_r_R_het_defined,
    n_r_R_het_undefined    = l5$n_r_R_het_undefined,
    level5_k2_degenerate   = l5$is_k2_degenerate,

    stringsAsFactors = FALSE
  )
}


# ============================================================================
# 3. Main Simulation Runner
# ============================================================================

run_simulation <- function(n_reps   = 100,
                            K_values = c(2, 3, 5, 10),
                            n_values = c(100, 200, 500, 1000),
                            p_values = c(2, 5, 10, 50),
                            scenarios = c("heteroscedastic", "spread_only",
                                          "homoscedastic", "nonlinear_var"),
                            verbose  = TRUE) {

  p_values <- sort(unique(p_values))
  p_max <- max(p_values)

  generators <- list(
    heteroscedastic = generate_heteroscedastic,
    spread_only     = generate_spread_only,
    homoscedastic   = generate_homoscedastic,
    nonlinear_var   = generate_nonlinear_variance
  )

  # Count total configurations
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

        # Generate one max-p template per replication, then subset it for
        # smaller p. This keeps the per-variable signal structure matched
        # across p-sensitivity comparisons.
        for (rep in 1:n_reps) {
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

  do.call(rbind, results[1:ri])
}


# ============================================================================
# 4. Run Simulations
# ============================================================================

cat("================================================================\n")
cat("I-WABA Simulation Study v3 — Matched p-Sensitivity + Level V Patch\n")
cat("================================================================\n\n")

sim_results <- run_simulation(
  n_reps   = 100,
  K_values = c(2, 3, 5, 10),
  n_values = c(100, 200, 500, 1000),
  p_values = c(2, 5, 10, 50),
  scenarios = c("heteroscedastic", "spread_only", "homoscedastic", "nonlinear_var"),
  verbose  = TRUE
)

saveRDS(sim_results, file = RESULTS_RDS)
write.csv(sim_results, file = RESULTS_CSV, row.names = FALSE)
cat(sprintf("\nSaved %d rows to %s / %s\n", nrow(sim_results), RESULTS_RDS, RESULTS_CSV))


# ============================================================================
# 5. Summary Tables
# ============================================================================

cat("\n================================================================\n")
cat("TABLE: Focal case n=500, p=10\n")
cat("================================================================\n")

focal <- sim_results %>% filter(n == 500, p == 10)

summary_table <- focal %>%
  group_by(scenario, K) %>%
  summarize(
    G_mean        = mean(info_gain),
    G_sd          = sd(info_gain),
    radius_prop   = mean(radius_prop),
    eta_B_waba    = mean(mean_eta_B_waba),
    eta_B_int     = mean(mean_eta_B_int),
    delta_eta     = mean(mean_delta_eta),
    n_reclass     = mean(n_reclassified),
    pi_R          = mean(mean_pi_R),
    E_R_het       = mean(mean_E_R_het_int),
    r_R_billard   = mean_or_na(mean_r_R_billard),
    r_R_het_signed = mean_or_na(mean_r_R_het_signed),
    r_R_het_pos   = mean_or_na(prop_r_R_het_positive),
    r_R_het_neg   = mean_or_na(prop_r_R_het_negative),
    r_R_het_undef = mean_or_na(prop_r_R_het_undefined),
    level5_k2_degenerate = any(level5_k2_degenerate),
    bw_int        = mean(mean_abs_bw_int),
    bw_waba       = mean(mean_abs_bw_waba),
    cons_int      = mean(consistency_rate_int),
    cons_waba     = mean(consistency_rate_waba),
    .groups = "drop"
  ) %>%
  mutate(
    level5_note = mapply(format_level5_note, K, r_R_het_pos, r_R_het_neg, r_R_het_undef),
    r_R_het_plot = if_else(K == 2, NA_real_, r_R_het_signed)
  )

# Level I
cat("\n--- Level I: Entity Analysis ---\n")
cat(sprintf("%-18s %3s %8s %7s %8s %9s %8s %6s %8s\n",
            "Scenario", "K", "G(mean)", "G(sd)", "RadProp",
            "eta_WABA", "eta_INT", "D.eta", "Reclass"))
for (i in 1:nrow(summary_table)) {
  r <- summary_table[i, ]
  cat(sprintf("%-18s %3d %8.3f %7.3f %8.3f %9.3f %8.3f %6.3f %8.1f\n",
              r$scenario, r$K, r$G_mean, r$G_sd, r$radius_prop,
              r$eta_B_waba, r$eta_B_int, r$delta_eta, r$n_reclass))
}

# Level IV
cat("\n--- Level IV: Dispersion-Source Diagnostics ---\n")
cat(sprintf("%-18s %3s %6s %8s\n", "Scenario", "K", "pi_R", "E_R,het"))
for (i in 1:nrow(summary_table)) {
  r <- summary_table[i, ]
  cat(sprintf("%-18s %3d %6.3f %8.3f\n", r$scenario, r$K, r$pi_R, r$E_R_het))
}

# Level V
cat("\n--- Level V: Dispersion Association ---\n")
cat(sprintf("%-18s %3s %10s %12s %s\n", "Scenario", "K", "r_R^Bill", "r_R,het", "Note"))
for (i in 1:nrow(summary_table)) {
  r <- summary_table[i, ]
  cat(sprintf("%-18s %3d %10s %12s %s\n",
              r$scenario, r$K,
              format_num_or_dash(r$r_R_billard),
              format_num_or_dash(r$r_R_het_signed),
              r$level5_note))
}

# Level II/III
cat("\n--- Level II/III: Functional Relationship & Consistency ---\n")
cat(sprintf("%-18s %3s %9s %10s %9s %10s\n",
            "Scenario", "K", "|Bw_INT|", "|Bw_WABA|", "Cons_INT", "Cons_WABA"))
for (i in 1:nrow(summary_table)) {
  r <- summary_table[i, ]
  cat(sprintf("%-18s %3d %9.3f %10.3f %9.3f %10.3f\n",
              r$scenario, r$K, r$bw_int, r$bw_waba, r$cons_int, r$cons_waba))
}


# ============================================================================
# 6. Publication Figures
# ============================================================================

cat("\n================================================================\n")
cat("Generating publication figures...\n")
cat("================================================================\n")

# Color palette
scenario_colors <- c(
  "heteroscedastic" = "#1f77b4",
  "spread_only"     = "#ff7f0e",
  "homoscedastic"   = "#2ca02c",
  "nonlinear_var"   = "#d62728"
)
scenario_labels <- c(
  "heteroscedastic" = "Heteroscedastic",
  "spread_only"     = "Spread-Only",
  "homoscedastic"   = "Homoscedastic",
  "nonlinear_var"   = "Nonlinear-Var"
)

# --- Figure 1: Information Gain vs K (two panels) ---
fig1_data <- summary_table %>%
  mutate(scenario_label = scenario_labels[scenario])

# Panel (a): het, homo, nonlinear
p1a <- ggplot(fig1_data %>% filter(scenario != "spread_only"),
              aes(x = K, y = G_mean, color = scenario_label, shape = scenario_label)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(title = "(a) Heteroscedastic, Homoscedastic, Nonlinear-Var",
       x = "K", y = "Information Gain G", color = "Scenario", shape = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# Panel (b): spread-only
p1b <- ggplot(fig1_data %>% filter(scenario == "spread_only"),
              aes(x = K, y = G_mean)) +
  geom_line(linewidth = 1, color = "#ff7f0e") +
  geom_point(size = 3, color = "#ff7f0e") +
  labs(title = "(b) Spread-Only (separate scale)",
       x = "K", y = "Information Gain G") +
  theme_bw(base_size = 11)

fig1 <- grid.arrange(p1a, p1b, ncol = 2, widths = c(1.3, 1))
ggsave(FIG_INFO_GAIN_FILE, fig1, width = 12, height = 5)

# --- Figure 2: Eta Improvement ---
summary_by_n <- focal %>%
  group_by(scenario, K, n) %>%
  summarize(delta_eta = mean(mean_delta_eta), .groups = "drop") %>%
  mutate(scenario_label = scenario_labels[scenario])

# Use the full grid for eta improvement
eta_grid <- sim_results %>%
  filter(p == 10, scenario %in% c("heteroscedastic", "spread_only")) %>%
  group_by(scenario, K, n) %>%
  summarize(delta_eta = mean(mean_delta_eta), .groups = "drop") %>%
  mutate(scenario_label = scenario_labels[scenario])

fig2 <- ggplot(eta_grid, aes(x = K, y = delta_eta, color = factor(n), shape = factor(n))) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  facet_wrap(~scenario_label, scales = "free_y") +
  labs(title = "Improvement in Mean Between-Group Eta",
       x = "K", y = expression(bar(Delta) * eta[B]),
       color = "n", shape = "n") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(FIG_ETA_IMPROVEMENT_FILE, fig2, width = 12, height = 5)

# --- Figure 3: Dispersion Diagnostics (two panels) ---
p3a <- ggplot(fig1_data, aes(x = K, y = pi_R, color = scenario_label, shape = scenario_label)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  labs(title = expression("(a) Mean Radius Share " * bar(pi)[R]),
       x = "K", y = expression(bar(pi)[R]),
       color = "Scenario", shape = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p3b <- ggplot(fig1_data, aes(x = K, y = E_R_het, color = scenario_label, shape = scenario_label)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  labs(title = expression("(b) Heterogeneity Index " * bar(E)[R*","*het]^"(Int)"),
       x = "K", y = expression(bar(E)[R*","*het]^"(Int)"),
       color = "Scenario", shape = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

fig3 <- grid.arrange(p3a, p3b, ncol = 2)
ggsave(FIG_DISPERSION_FILE, fig3, width = 12, height = 5)

# --- Figure 4: Heatmap (info gain, heteroscedastic, K vs p) ---
heatmap_data <- sim_results %>%
  filter(scenario == "heteroscedastic") %>%
  group_by(K, p) %>%
  summarize(G = mean(info_gain), .groups = "drop")

fig4 <- ggplot(heatmap_data, aes(x = factor(K), y = factor(p), fill = G)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", G)), size = 4) +
  scale_fill_gradient2(low = "white", mid = "steelblue", high = "darkblue",
                        midpoint = median(heatmap_data$G)) +
  labs(title = "Information Gain (Heteroscedastic)",
       x = "K", y = "p", fill = "G") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")
ggsave(FIG_HEATMAP_FILE, fig4, width = 8, height = 6)

# --- Figure 5: Radius Association (Level V) ---
fig5 <- ggplot(fig1_data %>% filter(K > 2),
               aes(x = K, y = r_R_het_plot, color = scenario_label, shape = scenario_label)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  labs(title = expression("Level V: Signed Heterogeneity-Based Radius Correlation " * bar(r)[R*","*het]),
       subtitle = "K = 2 omitted: sign-only / structurally degenerate diagnostic",
       x = "K", y = expression(bar(r)[R*","*het]),
       color = "Scenario", shape = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(FIG_RADIUS_FILE, fig5, width = 8, height = 5)


# ============================================================================
# 7. LaTeX Tables
# ============================================================================

cat("\n================================================================\n")
cat("LaTeX Tables\n")
cat("================================================================\n\n")

# Table 3: Level I
cat("% Table 3: Level I Entity Analysis\n")
cat("\\begin{table}[ht]\n\\centering\n")
cat("\\caption{Simulation results: Level~I entity analysis ($n=500$, $p=10$, 100 replications).}\n")
cat("\\label{tab:sim_results}\n")
cat("\\resizebox{\\textwidth}{!}{%\n")
cat("\\begin{tabular}{llcccccc}\n\\toprule\n")
cat("Scenario & $K$ & $\\mathcal{G}$ (Mean) & $\\mathcal{G}$ (SD) & $\\bar{\\pi}_R$ & ",
    "$\\bar{\\eta}_B^{\\text{WABA}}$ & $\\bar{\\eta}_B^{\\text{I-WABA}}$ & Reclass. \\\\\n",
    sep = "")
cat("\\midrule\n")

prev_sc <- ""
for (i in 1:nrow(summary_table)) {
  r <- summary_table[i, ]
  sc <- scenario_labels[r$scenario]
  if (prev_sc != "" && sc != prev_sc) cat("\\midrule\n")
  cat(sprintf("%s & %d & %.3f & %.3f & %.3f & %.3f & %.3f & %.1f \\\\\n",
              sc, r$K, r$G_mean, r$G_sd, r$radius_prop,
              r$eta_B_waba, r$eta_B_int, r$n_reclass))
  prev_sc <- sc
}
cat("\\bottomrule\n\\end{tabular}%\n}\n\\end{table}\n\n")

# Table: Level IV
cat("% Level IV: Dispersion-Source Diagnostics\n")
cat("\\begin{table}[ht]\n\\centering\n")
cat("\\caption{Level~IV dispersion-source diagnostics ($n=500$, $p=10$).}\n")
cat("\\label{tab:sim_level_iv}\n")
cat("\\begin{tabular}{llcc}\n\\toprule\n")
cat("Scenario & $K$ & $\\bar{\\pi}_R$ & $\\bar{E}_{R,\\text{het}}^{(\\text{Int})}$ \\\\\n")
cat("\\midrule\n")

prev_sc <- ""
for (i in 1:nrow(summary_table)) {
  r <- summary_table[i, ]
  sc <- scenario_labels[r$scenario]
  if (prev_sc != "" && sc != prev_sc) cat("\\midrule\n")
  cat(sprintf("%s & %d & %.3f & %.3f \\\\\n", sc, r$K, r$pi_R, r$E_R_het))
  prev_sc <- sc
}
cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")

# Table: Level V
cat("% Level V: Dispersion Association\n")
cat("\\begin{table}[ht]\n\\centering\n")
cat("\\caption{Level~V dispersion association ($n=500$, $p=10$). $\\bar{r}_R^{(\\text{Billard})}$ is reported on its natural $[0,1]$ scale. For $K=2$, $r_{R,\\text{het}}$ is structurally degenerate and therefore reported via a sign-only note rather than a graded mean.}\n")
cat("\\label{tab:sim_level_v}\n")
cat("\\begin{tabular}{llccc}\n\\toprule\n")
cat("Scenario & $K$ & $\\bar{r}_R^{(\\text{Billard})}$ & $\\bar{r}_{R,\\text{het}}$ & Note \\\\\n")
cat("\\midrule\n")

prev_sc <- ""
for (i in 1:nrow(summary_table)) {
  r <- summary_table[i, ]
  sc <- scenario_labels[r$scenario]
  if (prev_sc != "" && sc != prev_sc) cat("\\midrule\n")
  level5_note <- if (nzchar(r$level5_note)) r$level5_note else "--"
  cat(sprintf("%s & %d & %s & %s & %s \\\\\n",
              sc, r$K,
              format_num_or_dash(r$r_R_billard),
              format_num_or_dash(r$r_R_het_signed),
              level5_note))
  prev_sc <- sc
}
cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

cat("\n================================================================\n")
cat("Simulation study complete.\n")
cat("================================================================\n")
