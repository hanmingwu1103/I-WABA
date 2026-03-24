# Real-Data Run Log v2

## Milestone: Preflight

- Files changed
  - `iwaba_inference_helpers.R`
  - `real_data_bh1996.R`
  - `real_data_14cancer.R`
  - `real_data_analysis.R`
  - `results/realdata_execution_map_v2.md`
  - `results/realdata_runlog_v2.md`

- Commands run
  - `Get-ChildItem results -Force | Select-Object Mode,Length,LastWriteTime,Name | Format-Table -AutoSize | Out-String -Width 260`
  - `Get-Content real_data_analysis.R | Out-String -Width 260`
  - `Get-Content real_data_bh1996.R | Out-String -Width 260`
  - `Get-Content real_data_14cancer.R | Out-String -Width 260`
  - `Get-Content iwaba_functions.R | Out-String -Width 260`
  - `Get-Content iwaba_inference_helpers.R | Out-String -Width 260`
  - `rg -n "Top-50 branch|All-gene branch|Done\\." real_data_14cancer.R`
  - `& 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe' -e "invisible(parse(file='iwaba_inference_helpers.R')); invisible(parse(file='real_data_bh1996.R')); invisible(parse(file='real_data_14cancer.R')); invisible(parse(file='real_data_analysis.R')); cat('syntax-ok\n')"`

- Outputs created
  - `results/realdata_execution_map_v2.md`
  - `results/realdata_runlog_v2.md`

- Seeds used
  - No RNG executed in preflight
  - Fixed rerun seed standardized for later milestones: `20260321`

- What remains
  - Milestone 1: run `bh1996` with checkpointed saves and metadata
  - Milestone 2: run `14-Cancer` top-50 stage with restricted bootstrap (`B=399`, `max-pairs=10`)
  - Milestone 3: run `14-Cancer` all-gene stage with the same tag so it can write the sensitivity comparison
  - If desired, add `Rscript.exe` to `PATH`; current shell requires the full executable path

## Milestone: bh1996 / Task A - Primary SD-radius rerun

- Files changed
  - `bh1996_primary_v3.csv`
  - `bh1996_primary_v3.rds`
  - `fig_bh1996_analysis_v3.pdf`
  - `bh1996_metadata_v1.csv`
  - `bh1996_metadata_v1.rds`

- Commands run
  - `real_data_bh1996.R --B=399 --alpha=0.05 --seed=20260321 --output-dir=results`
  - `Subtask A: original-scale SD-radius analysis with within-group bootstrap B=399`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_primary_v3.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_primary_v3.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/fig_bh1996_analysis_v3.pdf`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_metadata_v1.rds`

- Seeds used
  - `Primary bootstrap seed = 20260321`

- What remains
  - Task B: group-size summary
  - Task C: IQR-radius sensitivity
  - Task D: scale-comparability sensitivity
  - Task E: interpretive inference table

## Milestone: bh1996 / Task B - Group-size summary

- Files changed
  - `bh1996_group_sizes_v1.csv`
  - `bh1996_group_sizes_v1.rds`
  - `bh1996_metadata_v1.csv`
  - `bh1996_metadata_v1.rds`

- Commands run
  - `real_data_bh1996.R --B=399 --alpha=0.05 --seed=20260321 --output-dir=results`
  - `Subtask B: company-size summary from table(GRP)`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_group_sizes_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_group_sizes_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_metadata_v1.rds`

- Seeds used
  - `No RNG used`

- What remains
  - Task C: IQR-radius sensitivity
  - Task D: scale-comparability sensitivity
  - Task E: interpretive inference table

## Milestone: bh1996 / Task C - IQR-radius sensitivity

- Files changed
  - `bh1996_radius_sensitivity_v1.csv`
  - `bh1996_radius_sensitivity_v1.rds`
  - `fig_bh1996_radius_sensitivity_v1.pdf`
  - `bh1996_metadata_v1.csv`
  - `bh1996_metadata_v1.rds`

- Commands run
  - `real_data_bh1996.R --B=399 --alpha=0.05 --seed=20260321 --output-dir=results`
  - `Subtask C: original-scale IQR-radius analysis with within-group bootstrap B=399`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_radius_sensitivity_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_radius_sensitivity_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/fig_bh1996_radius_sensitivity_v1.pdf`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_metadata_v1.rds`

- Seeds used
  - `IQR bootstrap seed = 20260322`

- What remains
  - Task D: scale-comparability sensitivity
  - Task E: interpretive inference table

## Milestone: bh1996 / Task D - Scale-comparability sensitivity

- Files changed
  - `bh1996_scale_sensitivity_v1.csv`
  - `bh1996_scale_sensitivity_v1.rds`
  - `bh1996_metadata_v1.csv`
  - `bh1996_metadata_v1.rds`

- Commands run
  - `real_data_bh1996.R --B=399 --alpha=0.05 --seed=20260321 --output-dir=results`
  - `Subtask D: globally standardized SD-radius rerun without bootstrap`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_scale_sensitivity_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_scale_sensitivity_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_metadata_v1.rds`

- Seeds used
  - `No RNG used`

- What remains
  - Task E: interpretive inference table

## Milestone: bh1996 / Task E - Interpretive inference table

- Files changed
  - `bh1996_inference_interpretive_v1.csv`
  - `bh1996_inference_interpretive_v1.rds`
  - `bh1996_metadata_v1.csv`
  - `bh1996_metadata_v1.rds`

- Commands run
  - `real_data_bh1996.R --B=399 --alpha=0.05 --seed=20260321 --output-dir=results`
  - `Subtask E: merge Brown-Forsythe, Level IV, and primary bootstrap delta_E intervals`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_inference_interpretive_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_inference_interpretive_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/bh1996_metadata_v1.rds`

- Seeds used
  - `Primary bootstrap seed reused for interpretive delta_E intervals = 20260321`

- What remains
  - bh1996 milestone complete

## Milestone: 14-Cancer / Task A - Existing top-50 and all-gene branches

- Files changed
  - `cancer14_top50_v3.csv`
  - `cancer14_top50_v3.rds`
  - `cancer14_allgenes_v3.csv`
  - `cancer14_allgenes_v3.rds`
  - `fig_14cancer_top50_v3.pdf`
  - `cancer14_metadata_v1.csv`
  - `cancer14_metadata_v1.rds`

- Commands run
  - `real_data_14cancer.R --stage=full --B=399 --alpha=0.05 --seed=20260321 --max-pairs=10 --output-dir=results`
  - `Task A top-50: top 50 ANOVA F-ranked probes with restricted bootstrap B=399`
  - `Task A all-genes: all retained probes with Brown-Forsythe screening only`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_top50_v3.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_top50_v3.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_allgenes_v3.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_allgenes_v3.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/fig_14cancer_top50_v3.pdf`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.rds`

- Seeds used
  - `Top-50 bootstrap seed = 20260321`
  - `All-gene branch uses no RNG`

- What remains
  - Task B: top-25 branch
  - Task C: moderate-F 50 branch
  - Task D: selection-sensitivity comparison
  - Task E: all-gene candidate-pair count

## Milestone: 14-Cancer / Task B - Top-25 branch

- Files changed
  - `cancer14_top25_v1.csv`
  - `cancer14_top25_v1.rds`
  - `cancer14_metadata_v1.csv`
  - `cancer14_metadata_v1.rds`

- Commands run
  - `real_data_14cancer.R --stage=full --B=399 --alpha=0.05 --seed=20260321 --max-pairs=10 --output-dir=results`
  - `Task B: top 25 ANOVA F-ranked probes with restricted bootstrap B=399`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_top25_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_top25_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.rds`

- Seeds used
  - `Top-25 bootstrap seed = 20260322`

- What remains
  - Task C: moderate-F 50 branch
  - Task D: selection-sensitivity comparison
  - Task E: all-gene candidate-pair count

## Milestone: 14-Cancer / Task C - Moderate-F 50 branch

- Files changed
  - `cancer14_moderate50_v1.csv`
  - `cancer14_moderate50_v1.rds`
  - `cancer14_selected_genes_moderate50_v1.csv`
  - `cancer14_selected_genes_moderate50_v1.rds`
  - `cancer14_metadata_v1.csv`
  - `cancer14_metadata_v1.rds`

- Commands run
  - `real_data_14cancer.R --stage=full --B=399 --alpha=0.05 --seed=20260321 --max-pairs=10 --output-dir=results`
  - `Task C: sample 50 probes uniformly from middle 40% ANOVA F ranks (4819-11245) with selection seed 20260323 and restricted bootstrap B=399`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_moderate50_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_moderate50_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_selected_genes_moderate50_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_selected_genes_moderate50_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.rds`

- Seeds used
  - `Moderate-F gene-selection seed = 20260323`
  - `Moderate-F bootstrap seed = 20260324`

- What remains
  - Task D: selection-sensitivity comparison
  - Task E: all-gene candidate-pair count

## Milestone: 14-Cancer / Task A - Existing top-50 and all-gene branches

- Files changed
  - `cancer14_top50_v3.csv`
  - `cancer14_top50_v3.rds`
  - `cancer14_allgenes_v3.csv`
  - `cancer14_allgenes_v3.rds`
  - `fig_14cancer_top50_v3.pdf`
  - `cancer14_metadata_v1.csv`
  - `cancer14_metadata_v1.rds`

- Commands run
  - `real_data_14cancer.R --stage=full --B=399 --alpha=0.05 --seed=20260321 --max-pairs=10 --output-dir=results`
  - `Task A top-50: top 50 ANOVA F-ranked probes with restricted bootstrap B=399`
  - `Task A all-genes: all retained probes with Brown-Forsythe screening only`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_top50_v3.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_top50_v3.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_allgenes_v3.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_allgenes_v3.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/fig_14cancer_top50_v3.pdf`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.rds`

- Seeds used
  - `Top-50 bootstrap seed = 20260321`
  - `All-gene branch uses no RNG`

- What remains
  - Task B: top-25 branch
  - Task C: moderate-F 50 branch
  - Task D: selection-sensitivity comparison
  - Task E: all-gene candidate-pair count

## Milestone: 14-Cancer / Task B - Top-25 branch

- Files changed
  - `cancer14_top25_v1.csv`
  - `cancer14_top25_v1.rds`
  - `cancer14_metadata_v1.csv`
  - `cancer14_metadata_v1.rds`

- Commands run
  - `real_data_14cancer.R --stage=full --B=399 --alpha=0.05 --seed=20260321 --max-pairs=10 --output-dir=results`
  - `Task B: top 25 ANOVA F-ranked probes with restricted bootstrap B=399`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_top25_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_top25_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.rds`

- Seeds used
  - `Top-25 bootstrap seed = 20260322`

- What remains
  - Task C: moderate-F 50 branch
  - Task D: selection-sensitivity comparison
  - Task E: all-gene candidate-pair count

## Milestone: 14-Cancer / Task C - Moderate-F 50 branch

- Files changed
  - `cancer14_moderate50_v1.csv`
  - `cancer14_moderate50_v1.rds`
  - `cancer14_selected_genes_moderate50_v1.csv`
  - `cancer14_selected_genes_moderate50_v1.rds`
  - `cancer14_metadata_v1.csv`
  - `cancer14_metadata_v1.rds`

- Commands run
  - `real_data_14cancer.R --stage=full --B=399 --alpha=0.05 --seed=20260321 --max-pairs=10 --output-dir=results`
  - `Task C: sample 50 probes uniformly from middle 40% ANOVA F ranks (4819-11245) with selection seed 20260323 and restricted bootstrap B=399`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_moderate50_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_moderate50_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_selected_genes_moderate50_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_selected_genes_moderate50_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.rds`

- Seeds used
  - `Moderate-F gene-selection seed = 20260323`
  - `Moderate-F bootstrap seed = 20260324`

- What remains
  - Task D: selection-sensitivity comparison
  - Task E: all-gene candidate-pair count

## Milestone: 14-Cancer / Task D - Selection-sensitivity comparison

- Files changed
  - `cancer14_selection_sensitivity_v1.csv`
  - `cancer14_selection_sensitivity_v1.rds`
  - `fig_14cancer_selection_sensitivity_v1.pdf`
  - `cancer14_metadata_v1.csv`
  - `cancer14_metadata_v1.rds`

- Commands run
  - `real_data_14cancer.R --stage=full --B=399 --alpha=0.05 --seed=20260321 --max-pairs=10 --output-dir=results`
  - `Task D: compare top-50, top-25, moderate-F 50, and all-gene screening summaries`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_selection_sensitivity_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_selection_sensitivity_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/fig_14cancer_selection_sensitivity_v1.pdf`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.rds`

- Seeds used
  - `No new RNG used beyond saved branch outputs`

- What remains
  - Task E: all-gene candidate-pair count

## Milestone: 14-Cancer / Task E - All-gene pair-count scope

- Files changed
  - `cancer14_allgenes_paircount_v1.csv`
  - `cancer14_allgenes_paircount_v1.rds`
  - `cancer14_metadata_v1.csv`
  - `cancer14_metadata_v1.rds`

- Commands run
  - `real_data_14cancer.R --stage=full --B=399 --alpha=0.05 --seed=20260321 --max-pairs=10 --output-dir=results`
  - `Task E: count BH-positive heterogeneous genes and choose(m, 2) candidate pairs for the all-gene branch`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_allgenes_paircount_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_allgenes_paircount_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.rds`

- Seeds used
  - `No RNG used`

- What remains
  - 14-Cancer screening/selection milestone complete

## Milestone: 14-Cancer / Pairwise Task A - Top-50 exhaustive BH

- Files changed
  - `cancer14_pairwise_fdr_top50_v1.csv`
  - `cancer14_pairwise_fdr_top50_v1.rds`
  - `cancer14_pairwise_fdr_top50_summary_v1.csv`
  - `cancer14_pairwise_fdr_top50_summary_v1.rds`
  - `cancer14_pairwise_fdr_top50_checkpoint_v1.rds`
  - `cancer14_metadata_v1.csv`
  - `cancer14_metadata_v1.rds`

- Commands run
  - `real_data_14cancer.R --stage=pairwise_fdr --B=399 --alpha=0.05 --seed=20260321 --checkpoint-every=25 --output-dir=results`
  - `Task A: exhaustive within-group bootstrap over all 1,225 top-50 pairs with B=399 and checkpoints every 25 resamples`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_pairwise_fdr_top50_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_pairwise_fdr_top50_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_pairwise_fdr_top50_summary_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_pairwise_fdr_top50_summary_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_pairwise_fdr_top50_checkpoint_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.rds`

- Seeds used
  - `Top-50 pairwise bootstrap seed = 20260321`

- What remains
  - Task B: top-25 exhaustive pairwise BH
  - Task C: restricted all-gene illustration among top 50 E_R_het-ranked BH-significant genes

## Milestone: 14-Cancer / Pairwise Task B - Top-25 exhaustive BH

- Files changed
  - `cancer14_pairwise_fdr_top25_v1.csv`
  - `cancer14_pairwise_fdr_top25_v1.rds`
  - `cancer14_pairwise_fdr_top25_checkpoint_v1.rds`
  - `cancer14_metadata_v1.csv`
  - `cancer14_metadata_v1.rds`

- Commands run
  - `real_data_14cancer.R --stage=pairwise_fdr --B=399 --alpha=0.05 --seed=20260321 --checkpoint-every=25 --output-dir=results`
  - `Task B: exhaustive within-group bootstrap over all 300 top-25 pairs with B=399 and checkpoints every 25 resamples`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_pairwise_fdr_top25_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_pairwise_fdr_top25_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_pairwise_fdr_top25_checkpoint_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.rds`

- Seeds used
  - `Top-25 pairwise bootstrap seed = 20260322`

- What remains
  - Task C: restricted all-gene illustration among top 50 E_R_het-ranked BH-significant genes

## Milestone: 14-Cancer / Pairwise Task C - All-gene illustration

- Files changed
  - `cancer14_allgenes_pairwise_illustration_v1.csv`
  - `cancer14_allgenes_pairwise_illustration_v1.rds`
  - `cancer14_allgenes_pairwise_illustration_checkpoint_v1.rds`
  - `cancer14_metadata_v1.csv`
  - `cancer14_metadata_v1.rds`

- Commands run
  - `real_data_14cancer.R --stage=pairwise_fdr --B=399 --alpha=0.05 --seed=20260321 --checkpoint-every=25 --output-dir=results`
  - `Task C: pairwise bootstrap over the top 50 BH-significant all-gene probes ranked by E_R_het^(Int) with B=399 and checkpoints every 25 resamples`

- Outputs created
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_allgenes_pairwise_illustration_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_allgenes_pairwise_illustration_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_allgenes_pairwise_illustration_checkpoint_v1.rds`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.csv`
  - `D:/08-MyProjects-new/04-I-WABA/05-I-WABA_20260330/I-WABA_Rproject/results/cancer14_metadata_v1.rds`

- Seeds used
  - `All-gene illustration selection is deterministic from saved screening outputs`
  - `All-gene illustration bootstrap seed = 20260323`

- What remains
  - Pairwise multiplicity-control milestone complete
