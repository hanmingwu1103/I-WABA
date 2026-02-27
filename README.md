# I-WABA: Interval-based Within-And-Between Analysis

R implementation of I-WABA, an extension of classical WABA (Within-And-Between Analysis) that incorporates within-group spread information using Symbolic Data Analysis (SDA).

## Reference

Wu, H.-M. (2026). When Means Are Not Enough: Enriching Within-And-Between Analysis with Group Intervals. *Technical Report*.

## Method

Classical WABA decomposes total variance using group means only. I-WABA represents each group as an **interval-valued symbol** (center = group mean, radius = group SD) and applies SDA covariance to capture both location and spread differences between groups.

Key quantities:
- **Information Gain** (*G*): ratio of I-WABA between-group variance to classical between-group variance. *G* > 1 indicates I-WABA captures additional between-group information from heteroscedasticity.
- **Radius Proportion**: fraction of between-SDA variance attributable to group spread differences.
- **I-WABA eta** (*eta_B*, *eta_W*): correlation ratios computed under the SDA framework.

## Repository Structure

```
iwaba_functions.R        # Core I-WABA functions (Entity & Functional Analysis)
simulation_study.R       # Simulation study (4 scenarios x varying K, n, p)
real_data_analysis.R     # Real data analysis (Yammarino & Markham; GCM cancer)
real_data_14cancer.R     # 14-Cancer gene expression analysis (Ramaswamy et al.)
real_data_bh1996.R       # bh1996 multilevel dataset analysis (Bliese & Halverson)
data/                    # Datasets (14-cancer, bh1996)
```

## Usage

```r
source("iwaba_functions.R")

# Run I-WABA analysis
result <- iwaba_analysis(X, group_labels)
print_iwaba_summary(result)
```

## Dependencies

`MASS`, `ggplot2`, `reshape2`, `gridExtra`, `dplyr`, `xtable`

## Author

Han-Ming Wu, National Chengchi University
