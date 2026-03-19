# I-WABA: Interval-based Within-And-Between Analysis

R implementation of I-WABA, an extension of classical WABA (Within-And-Between Analysis) that incorporates within-group dispersion using Billard's endpoint covariance for interval-valued data.

## Reference

Wu, H.-M. (2026). I-WABA: Interval-Based Within-And-Between Analysis for Heterogeneous Group Dispersion. *Statistical Analysis and Data Mining: The ASA Data Science Journal*.

## Method

Classical WABA decomposes total variance using group means only. I-WABA represents each group as an **interval-valued symbol** (center = group mean, radius = group SD) and applies Billard's endpoint covariance (Eq. 21) to capture both location and spread differences between groups.

### Five-Level Analysis Framework

- **Level I (Entity Analysis):** Eta correlation ratios, E-test classification, information gain *G*, radius proportion.
- **Level II (Functional Relationship):** Weighted correlation decomposition using augmented total, between-group interval, and within-group correlations.
- **Level III (Consistency Check):** Agreement between Level I entity classification and Level II relationship dominance.
- **Level IV (Dispersion-Source Diagnostics):** Center-only (*E_C*) vs radius-only (*E_R*) attribution, heterogeneity-focused companion (*E_{R,het}*) using centered radii.
- **Level V (Dispersion Association):** Billard radius association (*r_R^Billard*) and heterogeneity-based radius correlation (*r_{R,het}*).

### Key Formula (Billard's Endpoint Covariance)

```
V_R^(j) = (1/3) * sum_k w_k * s_kj^2
```

This uses the **uncentered second moment** of group radii (Eq. 21), not the centered covariance `Cov(R,R')` (Eq. 15). The Level IV diagnostic `V_{R,het}` uses centered radii to isolate genuine heteroscedasticity from magnitude effects.

## Repository Structure

```
iwaba_functions.R        # Core I-WABA functions (Level I–V)
simulation_study.R       # Simulation study (4 scenarios x K x n x p, Level I–V)
real_data_bh1996.R       # bh1996 military dataset (K=99, Level I–V)
real_data_14cancer.R     # 14-Cancer gene expression (K=14, Level I–V)
real_data_analysis.R     # (Deprecated — see individual scripts above)
data/                    # Datasets (14-cancer, bh1996)
```

## Usage

```r
source("iwaba_functions.R")

# Full Level I–V analysis
result <- iwaba_full(X, group_labels)
print_iwaba_summary(result)

# Access individual levels
result$level1  # Entity analysis (info gain, eta, E-test)
result$level2  # Functional relationships (weighted correlations)
result$level3  # Consistency check
result$level4  # Dispersion-source diagnostics (pi_R, E_R_het)
result$level5  # Dispersion association (r_R_Billard, r_R_het)

# Legacy wrapper (backward compatible)
result <- iwaba_analysis(X, group_labels)
```

## Dependencies

`MASS`, `ggplot2`, `dplyr`, `tidyr`, `gridExtra`

## Author

Han-Ming Wu, National Chengchi University
