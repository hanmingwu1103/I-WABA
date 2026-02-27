################################################################################
## I-WABA Real Data Analysis: Bliese-Halverson (bh1996) Dataset
## Author: Han-Ming Wu (National Chengchi University)
##
## Data: bh1996 from the 'multilevel' R package
##   n = 7,382 soldiers, K = 99 army companies, p = 4 variables
##   Reference: Bliese & Halverson (1996), The Leadership Quarterly
################################################################################

source("R/iwaba_functions.R")

# ============================================================================
# 1. Load Data
# ============================================================================
cat("Loading bh1996 data...\n")

# Option A: From multilevel package
if (require(multilevel, quietly = TRUE)) {
  data(bh1996)
  df <- bh1996
  cat("  Loaded from multilevel R package\n")
} else {
  # Option B: From CSV export
  df <- read.csv("data/bh1996_rdata.csv", row.names = 1)
  cat("  Loaded from data/bh1996_rdata.csv\n")
}

variables <- c("COHES", "LEAD", "WBEING", "HRS")
group_col <- "GRP"

cat(sprintf("  n = %d, K = %d, p = %d\n",
            nrow(df), length(unique(df[[group_col]])), length(variables)))

# ============================================================================
# 2. I-WABA I: Entity Analysis
# ============================================================================
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("I-WABA I: ENTITY ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════\n")

X <- as.matrix(df[, variables])
g <- df[[group_col]]

result <- iwaba_analysis(X, g)
print_iwaba_summary(result)

# Print detailed table
cat("\nDetailed eta comparison:\n")
cat(sprintf("  %-8s %8s %8s %8s %8s %8s %8s %8s %8s\n",
            "Variable", "η_B(W)", "η_W(W)", "E(W)",
            "η_B(I)", "η_W(I)", "E(I)", "Δη_B", "E_ratio"))
for (i in 1:nrow(result$eta_results)) {
  r <- result$eta_results[i, ]
  cat(sprintf("  %-8s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.3f\n",
              r$variable, r$eta_B_waba, r$eta_W_waba, r$E_waba,
              r$eta_B_iwaba, r$eta_W_iwaba, r$E_iwaba,
              r$delta_eta, r$E_ratio))
}

# ============================================================================
# 3. I-WABA II: Functional Relationship Analysis
# ============================================================================
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("I-WABA II: FUNCTIONAL RELATIONSHIP ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Between-group correlations (from group means)
grp_means <- aggregate(. ~ GRP, data = df[, c(group_col, variables)], FUN = mean)
between_corr <- cor(grp_means[, variables])

# Within-group correlations (from residuals)
residuals <- df[, variables]
for (v in variables) {
  grp_mean <- ave(df[[v]], df[[group_col]], FUN = mean)
  residuals[[v]] <- df[[v]] - grp_mean
}
within_corr <- cor(residuals)

# SDA between-group correlations
sda_std <- sqrt(diag(result$sda$cov_sda))
sda_corr <- result$sda$cov_sda / outer(sda_std, sda_std)

cat("\nPairwise WABA II / I-WABA II:\n")
cat(sprintf("  %-15s %8s %8s %8s %10s %10s %10s %10s\n",
            "Pair", "r_B", "r_W", "r_SDA_B",
            "W2_B", "W2_W", "IW2_B", "IW2_W"))

pairs <- combn(seq_along(variables), 2)
for (col in 1:ncol(pairs)) {
  i <- pairs[1, col]
  j <- pairs[2, col]
  v1 <- variables[i]; v2 <- variables[j]

  r_b <- between_corr[v1, v2]
  r_w <- within_corr[v1, v2]
  r_sda <- sda_corr[i, j]

  # WABA II weighted components
  w2_b <- result$eta_results$eta_B_waba[i] * result$eta_results$eta_B_waba[j] * r_b
  w2_w <- result$eta_results$eta_W_waba[i] * result$eta_results$eta_W_waba[j] * r_w

  # I-WABA II weighted components
  iw2_b <- result$eta_results$eta_B_iwaba[i] * result$eta_results$eta_B_iwaba[j] * r_sda
  iw2_w <- result$eta_results$eta_W_iwaba[i] * result$eta_results$eta_W_iwaba[j] * r_w

  cat(sprintf("  %-15s %8.3f %8.3f %8.3f %10.4f %10.4f %10.4f %10.4f\n",
              paste0(v1, "--", v2), r_b, r_w, r_sda, w2_b, w2_w, iw2_b, iw2_w))
}

cat("\nDone!\n")
