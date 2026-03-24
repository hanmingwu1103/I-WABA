# Real-Data Execution Map v2

## Entry points

- `bh1996`: [`real_data_bh1996.R`](../real_data_bh1996.R)
  - Loads `bh1996` from `multilevel` or `data/bh1996_rdata.csv`
  - Runs `iwaba_full(X, g)` for Levels I-V
  - Runs `iwaba_brown_forsythe_matrix()` plus `iwaba_within_group_bootstrap_vector()` via `bh_boot_stat_fn`

- `14-Cancer`: [`real_data_14cancer.R`](../real_data_14cancer.R)
  - Loads `data/14cancer.*`
  - Preprocesses with variance filter, z-score scaling, and ANOVA `F_scores`
  - `--stage=top50` runs `iwaba_full(X_top50, cancer_labels)` plus restricted bootstrap inference
  - `--stage=allgenes` runs `iwaba_level1()` and `iwaba_level4()` on all retained probes plus Brown-Forsythe screening

- Runner guide: [`real_data_analysis.R`](../real_data_analysis.R)
  - Documents the staged commands

## Save paths

- Versioned outputs now save under `results/` through shared helpers in [`iwaba_inference_helpers.R`](../iwaba_inference_helpers.R):
  - `iwaba_save_table_outputs()` writes paired `.csv` and `.rds`
  - `iwaba_open_pdf_device()` writes `.pdf`
  - `iwaba_save_metadata()` writes run metadata as `.csv` and `.rds`
- Optional `--tag=` appends a run-specific suffix to each versioned filename.
- Historical outputs already present in `results/` use the older untagged naming scheme.

## Control knobs

- Radius type
  - Controlled in [`iwaba_functions.R`](../iwaba_functions.R)
  - `iwaba_level1()` defines `V_R = (1/3) * sum_k w_k * s_kj^2` using Billard endpoint covariance
  - `iwaba_level4()` defines `pi_R`, `V_R_het`, and `E_R_het_int`

- Bootstrap resamples
  - `real_data_bh1996.R`: `--B=` -> `BOOT_B`
  - `real_data_14cancer.R`: `--B=` -> `BOOT_B`
  - Both pass `BOOT_B` into `iwaba_within_group_bootstrap_vector()`

- Feature selection
  - `bh1996`: none
  - `14-Cancer`: `gene_var > 1e-10`, `X_scaled <- scale(X_all)`, `F_scores`, `top50_idx`

- Pairwise inference
  - `bh1996`: `bh_boot_stat_fn` bootstraps all Level II / Level V pairs
  - `14-Cancer top50`: `MAX_PAIR_BOOT` and `selected_pair_labels` restrict pairwise bootstrap scope
  - `14-Cancer allgenes`: no all-gene pairwise bootstrap; the stage can read the saved top-50 summary to write `cancer_sensitivity_comparison_*`

## Staged commands

1. `& 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe' real_data_bh1996.R --B=399 --alpha=0.05 --seed=20260321 --output-dir=results --tag=realdata_v2`
2. `& 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe' real_data_14cancer.R --stage=top50 --B=399 --alpha=0.05 --seed=20260321 --max-pairs=10 --output-dir=results --tag=realdata_v2`
3. `& 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe' real_data_14cancer.R --stage=allgenes --alpha=0.05 --seed=20260321 --output-dir=results --tag=realdata_v2`

Use the same `--tag=realdata_v2` for milestones 2 and 3 so the all-gene stage can reuse the saved top-50 summary when it writes the sensitivity comparison.
