################################################################################
## I-WABA Real Data Analysis: 14-Cancer Gene Expression Dataset
## Author: Han-Ming Wu (National Chengchi University)
##
## Data: Ramaswamy et al. (2001) 14-Cancer benchmark
##   n = 198 tumor samples, K = 14 cancer types, p = 16,063 genes
##   Source: https://hastie.su.domains/ElemStatLearn/datasets/
##   Reference: Ramaswamy et al. (2001), PNAS 98(26): 15149-15154
################################################################################

source("R/iwaba_functions.R")

CANCER_NAMES <- c(
  "Breast", "Prostate", "Lung", "Colorectal", "Lymphoma",
  "Bladder", "Melanoma", "Uterus", "Leukemia", "Renal",
  "Pancreas", "Ovary", "Mesothelioma", "CNS"
)

# ============================================================================
# 1. Load and Merge Train/Test Data
# ============================================================================
cat("═══════════════════════════════════════════════════════════════\n")
cat("Loading 14-Cancer data (Ramaswamy et al., 2001, PNAS)\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Read space-delimited files (genes x samples)
xtrain <- as.matrix(read.table("data/14cancer.xtrain"))
xtest  <- as.matrix(read.table("data/14cancer.xtest"))

ytrain <- scan("data/14cancer.ytrain", what = integer(), quiet = TRUE)
ytest  <- scan("data/14cancer.ytest", what = integer(), quiet = TRUE)

cat(sprintf("  Training: %d genes x %d samples\n", nrow(xtrain), ncol(xtrain)))
cat(sprintf("  Test:     %d genes x %d samples\n", nrow(xtest), ncol(xtest)))

# Merge and transpose (samples x genes)
X_all <- t(cbind(xtrain, xtest))
y_all <- c(ytrain, ytest)
cancer_labels <- CANCER_NAMES[y_all]

n <- nrow(X_all)
p_full <- ncol(X_all)
K <- length(unique(y_all))

cat(sprintf("  Combined: n=%d, p=%d, K=%d\n", n, p_full, K))
cat(sprintf("  p/n ratio = %.1f (HDLSS)\n", p_full / n))

cat("\n  Cancer type distribution:\n")
for (i in 1:14) {
  cnt <- sum(y_all == i)
  cat(sprintf("    %2d. %-15s: n=%d\n", i, CANCER_NAMES[i], cnt))
}

# ============================================================================
# 2. Preprocessing: Z-score + ANOVA F-score Top-50
# ============================================================================
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("Preprocessing: Z-score normalization + ANOVA top-50 genes\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Remove zero-variance genes
gene_var <- apply(X_all, 2, var)
X_all <- X_all[, gene_var > 1e-10]
cat(sprintf("  After QC: %d genes\n", ncol(X_all)))

# Z-score normalization (gene-wise)
X_scaled <- scale(X_all)

# ANOVA F-score for each gene
cat("  Computing ANOVA F-scores...\n")
grand_mean <- colMeans(X_scaled)
BSS <- numeric(ncol(X_scaled))
WSS <- numeric(ncol(X_scaled))

for (ct in unique(y_all)) {
  mask <- y_all == ct
  n_k <- sum(mask)
  group_mean <- colMeans(X_scaled[mask, , drop = FALSE])
  BSS <- BSS + n_k * (group_mean - grand_mean)^2
  WSS <- WSS + colSums(sweep(X_scaled[mask, , drop = FALSE], 2, group_mean)^2)
}

df_between <- K - 1
df_within <- n - K
F_scores <- (BSS / df_between) / (WSS / df_within)
F_scores[is.nan(F_scores)] <- 0

# Select top 50
top50_idx <- order(F_scores, decreasing = TRUE)[1:50]
X_top50 <- X_scaled[, top50_idx]
colnames(X_top50) <- paste0("Gene_", 1:50)

cat(sprintf("  Top-50 F-score range: %.2f to %.2f\n",
            min(F_scores[top50_idx]), max(F_scores[top50_idx])))
cat(sprintf("  Final matrix: %d x %d\n", nrow(X_top50), ncol(X_top50)))

# ============================================================================
# 3. I-WABA Analysis
# ============================================================================
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("I-WABA Analysis\n")
cat("═══════════════════════════════════════════════════════════════\n")

result <- iwaba_analysis(X_top50, cancer_labels)
print_iwaba_summary(result)

# Classification
result$eta_results$class_waba <- classify_e(result$eta_results$E_waba)
result$eta_results$class_iwaba <- classify_e(result$eta_results$E_iwaba)

reclassified <- sum(result$eta_results$class_waba != result$eta_results$class_iwaba)
cat(sprintf("\n  Reclassified genes: %d / 50\n", reclassified))

cat("\n  WABA I classification:\n")
print(table(result$eta_results$class_waba))
cat("\n  I-WABA I classification:\n")
print(table(result$eta_results$class_iwaba))

# ============================================================================
# 4. Heteroscedasticity Analysis
# ============================================================================
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("Heteroscedasticity Analysis\n")
cat("═══════════════════════════════════════════════════════════════\n")

mean_radius <- rowMeans(result$sda$radii)
names(mean_radius) <- result$sda$group_names
cat("\n  Mean within-group SD by cancer type:\n")
for (ct in names(sort(mean_radius, decreasing = TRUE))) {
  cat(sprintf("    %-15s: %.4f\n", ct, mean_radius[ct]))
}

# ============================================================================
# 5. Generate Figures
# ============================================================================
cat("\n  Generating figures...\n")

pdf("R/fig_14cancer_analysis.pdf", width = 14, height = 12)
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))

# (a) Eta scatter plot
plot(result$eta_results$eta_B_waba, result$eta_results$eta_B_iwaba,
     xlab = expression(eta[B]^"(WABA)"),
     ylab = expression(eta[B]^"(I-WABA)"),
     main = "(a) Entity Analysis: I-WABA vs WABA",
     pch = 19, col = "steelblue", cex = 1.2)
abline(0, 1, col = "red", lty = 2, lwd = 2)
grid()

# (b) Delta eta histogram
hist(result$eta_results$delta_eta, breaks = 20, col = "#DD8452",
     border = "black", main = expression("(b) Distribution of " * Delta * eta[B]),
     xlab = expression(Delta * eta[B]))
abline(v = mean(result$eta_results$delta_eta), col = "blue", lwd = 2)
abline(v = 0, col = "red", lty = 2)

# (c) Heteroscedasticity boxplot
sorted_radii <- result$sda$radii[order(mean_radius, decreasing = TRUE), ]
boxplot(t(sorted_radii), las = 2, col = "#55A868",
        main = "(c) Heteroscedasticity by Cancer Type",
        ylab = "Within-group SD", cex.axis = 0.7)

# (d) E-ratio scatter
plot(result$eta_results$E_waba, result$eta_results$E_iwaba,
     xlab = expression(E^"(WABA)"), ylab = expression(E^"(I-WABA)"),
     main = "(d) E-test Ratios: I-WABA vs WABA",
     pch = 19, col = "steelblue", cex = 1.2)
abline(0, 1, col = "red", lty = 2, lwd = 2)
abline(h = 1, v = 1, col = "gray", lty = 3)
grid()

dev.off()
cat("  Saved: R/fig_14cancer_analysis.pdf\n")

# ============================================================================
# 6. Sensitivity Analysis: All Genes (p = 16,063)
# ============================================================================
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("Sensitivity Analysis: I-WABA on ALL genes (p = 16,063)\n")
cat("═══════════════════════════════════════════════════════════════\n")

result_all <- iwaba_analysis(X_scaled, cancer_labels)
print_iwaba_summary(result_all)

# Classification
result_all$eta_results$class_waba <- classify_e(result_all$eta_results$E_waba)
result_all$eta_results$class_iwaba <- classify_e(result_all$eta_results$E_iwaba)

reclass_all <- sum(result_all$eta_results$class_waba != result_all$eta_results$class_iwaba)
p_all <- ncol(X_scaled)
cat(sprintf("\n  Reclassified genes: %d / %d (%.1f%%)\n",
            reclass_all, p_all, 100 * reclass_all / p_all))

cat(sprintf("  Proportion delta_eta > 0: %.4f (%d/%d)\n",
            mean(result_all$eta_results$delta_eta > 0),
            sum(result_all$eta_results$delta_eta > 0), p_all))

cat("\n  WABA I classification (all genes):\n")
print(table(result_all$eta_results$class_waba))
cat("\n  I-WABA I classification (all genes):\n")
print(table(result_all$eta_results$class_iwaba))

# Reclassification transitions
if (reclass_all > 0) {
  changed <- result_all$eta_results$class_waba != result_all$eta_results$class_iwaba
  transitions <- table(
    From = result_all$eta_results$class_waba[changed],
    To = result_all$eta_results$class_iwaba[changed]
  )
  cat("\n  Reclassification transitions:\n")
  print(transitions)
}

# Comparison summary
cat("\n  ══════════════════════════════════════════════════════════\n")
cat("  TOP-50 vs ALL-GENE COMPARISON\n")
cat("  ══════════════════════════════════════════════════════════\n")
cat(sprintf("  %-30s %10s %12s\n", "Metric", "Top-50", "All genes"))
cat(sprintf("  %-30s %10d %12d\n", "Number of genes", 50, p_all))
cat(sprintf("  %-30s %10.4f %12.4f\n", "Information Gain G",
            result$info_gain, result_all$info_gain))
cat(sprintf("  %-30s %10.4f %12.4f\n", "Radius Proportion",
            result$radius_prop, result_all$radius_prop))
cat(sprintf("  %-30s %10.4f %12.4f\n", "Mean eta_B (WABA)",
            mean(result$eta_results$eta_B_waba),
            mean(result_all$eta_results$eta_B_waba)))
cat(sprintf("  %-30s %10.4f %12.4f\n", "Mean eta_B (I-WABA)",
            mean(result$eta_results$eta_B_iwaba),
            mean(result_all$eta_results$eta_B_iwaba)))
cat(sprintf("  %-30s %10.4f %12.4f\n", "Mean delta eta_B",
            mean(result$eta_results$delta_eta),
            mean(result_all$eta_results$delta_eta)))

cat("\nDone!\n")
