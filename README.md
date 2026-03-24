# I-WABA: Interval-Based Within-And-Between Analysis

This repository contains the reproducible R workflows and manuscript materials
for the I-WABA project. The current snapshot matches the fourth-revision (`R4`)
manuscript bundle and the archived workflow used for the latest JRSSB revision
cycle.

## Repository Contents

- `iwaba_functions.R`: core I-WABA functions and summary helpers.
- `iwaba_inference_helpers.R`: bootstrap and inference utilities.
- `simulation_study.R`: main decomposition simulation study.
- `simulation_inference.R`: fixed-`K` inferential simulation study.
- `simulation_inference_robustness.R`: robustness study under skewness and imbalance.
- `simulation_correlation_ordering_summary.R`: summary script for the covariance-vs-correlation ordering discussion.
- `real_data_bh1996.R`: analysis of the `bh1996` multilevel benchmark dataset.
- `real_data_14cancer.R`: analysis of the 14-Cancer gene-expression dataset.
- `real_data_analysis.R`: legacy wrapper retained for convenience.
- `data/`: datasets used by the real-data analyses.
- `results/`: generated tables, summaries, figures, and supporting result objects.
- `LaTeX/`: curated manuscript sources and figure files for the current `R4` paper bundle.

## Reproducing the Main Results

Run the scripts below from the repository root as needed.

1. `source("simulation_study.R")`
2. `source("simulation_inference.R")`
3. `source("simulation_inference_robustness.R")`
4. `source("simulation_correlation_ordering_summary.R")`
5. `source("real_data_bh1996.R")`
6. `source("real_data_14cancer.R")`

The scripts write their main outputs to `results/`. The committed files in
`results/` are the manuscript versions used for this release snapshot.

## Dependencies

The main analyses rely on these R packages:

- `MASS`
- `ggplot2`
- `dplyr`
- `tidyr`
- `gridExtra`
- `multilevel` (tested with version `2.7.1`)

## Manuscript Materials

The current manuscript bundle is provided in `LaTeX/`, including:

- `LaTeX/I-WABA_Wu_R4.tex`
- `LaTeX/I-WABA_Wu_R4.bib`
- `LaTeX/I-WABA_Wu_R4.pdf`

The curated figure PDFs needed by `I-WABA_Wu_R4.tex` are included in the same
folder so the manuscript can be rebuilt from the release snapshot.

## Notes

- The archived JRSSB revision workflow used `multilevel` version `2.7.1`, matching the manuscript bibliography.
- The GitHub release history provides versioned snapshots of the code and manuscript materials used during the revision process.
