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
cat("  Milestone 1: Rscript real_data_bh1996.R --B=399 --alpha=0.05 --seed=20260321 --output-dir=results --tag=realdata_v2\n")
cat("  Outputs: results/bh1996_level1_v2_<tag>.csv/.rds, results/bh1996_level2_v2_<tag>.csv/.rds,\n")
cat("           results/bh1996_level3_v2_<tag>.csv/.rds, results/bh1996_level4_v2_<tag>.csv/.rds,\n")
cat("           results/bh1996_level5_v2_<tag>.csv/.rds, results/bh1996_summary_v2_<tag>.csv/.rds,\n")
cat("           results/bh1996_inference_v1_<tag>.csv/.rds, results/fig_bh1996_analysis_v2_<tag>.pdf,\n")
cat("           results/bh1996_metadata_v2_<tag>.csv/.rds\n\n")

cat("Dataset 2: 14-Cancer\n")
cat("  Milestone 2: Rscript real_data_14cancer.R --stage=top50 --B=399 --alpha=0.05 --seed=20260321 --max-pairs=10 --output-dir=results --tag=realdata_v2\n")
cat("  Milestone 3: Rscript real_data_14cancer.R --stage=allgenes --alpha=0.05 --seed=20260321 --output-dir=results --tag=realdata_v2\n")
cat("  Outputs: results/cancer_top50_level1_v2_<tag>.csv/.rds, results/cancer_top50_level2_v2_<tag>.csv/.rds,\n")
cat("           results/cancer_top50_level2_aggregate_v2_<tag>.csv/.rds, results/cancer_top50_level3_v2_<tag>.csv/.rds,\n")
cat("           results/cancer_top50_level4_v2_<tag>.csv/.rds, results/cancer_top50_level5_v2_<tag>.csv/.rds,\n")
cat("           results/cancer_top50_summary_v2_<tag>.csv/.rds, results/cancer_top50_screening_v1_<tag>.csv/.rds,\n")
cat("           results/cancer_top50_inference_v1_<tag>.csv/.rds, results/cancer_allgenes_summary_v2_<tag>.csv/.rds,\n")
cat("           results/cancer_allgenes_screening_v1_<tag>.csv/.rds, results/cancer_sensitivity_comparison_v2_<tag>.csv/.rds,\n")
cat("           results/fig_14cancer_top50_v2_<tag>.pdf, results/cancer_top50_metadata_v2_<tag>.csv/.rds,\n")
cat("           results/cancer_allgenes_metadata_v2_<tag>.csv/.rds\n")
