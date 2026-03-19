################################################################################
## I-WABA Real Data Analysis Entry Point
##
## This file is a lightweight runner guide for the revised real-data pipeline.
## Use the dataset-specific scripts directly so their versioned outputs remain
## explicit and reproducible.
################################################################################

cat("Revised real-data pipeline\n")
cat("==========================\n\n")

cat("Dataset 1: Bliese-Halverson (bh1996)\n")
cat("  Rscript real_data_bh1996.R\n")
cat("  Optional: Rscript real_data_bh1996.R --B=999 --alpha=0.05 --seed=20260319\n")
cat("  Outputs: bh1996_level1_v2.csv, bh1996_level2_v2.csv, bh1996_level3_v2.csv,\n")
cat("           bh1996_level4_v2.csv, bh1996_level5_v2.csv, bh1996_summary_v2.csv,\n")
cat("           bh1996_inference_v1.csv, fig_bh1996_analysis_v2.pdf\n\n")

cat("Dataset 2: 14-Cancer\n")
cat("  Rscript real_data_14cancer.R\n")
cat("  Optional: Rscript real_data_14cancer.R --B=999 --alpha=0.05 --seed=20260319 --max-pairs=10\n")
cat("  Outputs: cancer_top50_level1_v2.csv, cancer_top50_level2_v2.csv,\n")
cat("           cancer_top50_level2_aggregate_v2.csv, cancer_top50_level3_v2.csv,\n")
cat("           cancer_top50_level4_v2.csv, cancer_top50_level5_v2.csv,\n")
cat("           cancer_top50_summary_v2.csv, cancer_top50_screening_v1.csv,\n")
cat("           cancer_top50_inference_v1.csv, cancer_allgenes_summary_v2.csv,\n")
cat("           cancer_allgenes_screening_v1.csv, cancer_sensitivity_comparison_v2.csv,\n")
cat("           fig_14cancer_top50_v2.pdf\n")
