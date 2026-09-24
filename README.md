# Cancer sex disparities attributable to lifestyle factors and health conditions

This repository provides a concise R demonstration of the main calculations employed in the accompanying study: **Li ZY, Tan YT, Shen QM, Yang DN, Li HL, Xiang YB. Disparities in Cancer Incidence by Sex and Associated Lifestyle Factors and Health Conditions. ***JAMA Oncol***. Published online September 24, 2026. doi:10.1001/jamaoncol.2026.3747**. 
If you use these codes in your research, please cite our paper. 

## Repository contents
- male-to-female incidence rate ratios (IRRs) from Poisson regression;
- directly age-standardized incidence rates (ASIRs) using Segi-Doll weights;
- the male-female ASIR difference and 95% CIs;
- 5-year single-factor PAFs and joint PAFs from Cox models;
- attributable rate differences (ARDs);
- attributable proportions (APs); and
- percentile bootstrap confidence intervals for ARDs.

The original participant-level data cannot be shared because of data-protection and ethical restrictions. The repository therefore generates a mock dataset. The synthetic data are intended only to illustrate the calculations and do not reproduce the published results.

## Scripts Description
* `00_Generate_Mock_Data.R`: generates 10,000 simulated participants with both sexes, 13 cancer-incidence endpoints (pan-cancer plus 12 sites), covariates, and the risk factors used in the calculations. Mock event and exposure frequencies are set for a workable statistical demonstration, not to reproduce the original study's counts or precision.
* `01_IRR_ASIR.R`: calculates Poisson IRRs, direct ASIRs, and the male-female ASIR difference.
* `02_PAF_ARD_AP.R`: calculates sex-specific single-factor PAFs, joint PAFs, ARDs, and APs.
* `03_Bootstrap_ARD_CI.R`: calculates percentile bootstrap CIs for ARDs using 100 nonparametric bootstrap replicates.
* `AFcoxph_joint.R`: computes a joint PAF by setting all specified binary risk factors to their low-risk reference level simultaneously.

## Requirements
R (tested with R 4.5.3) and the following packages:

```r
install.packages(c("survival", "AF", "epitools"))
```

## Running the scripts
From the repository directory:

```r
source("00_Generate_Mock_Data.R")
source("01_IRR_ASIR.R")
source("02_PAF_ARD_AP.R")
source("03_Bootstrap_ARD_CI.R")
```

Run the scripts in this order. The bootstrap script reads the point estimates written by `02_PAF_ARD_AP.R` and does not recalculate them.

The scripts create:

```text
mock_data.rds
IRR_ASIR_results.csv
PAF_ARD_AP_results.csv
ARD_bootstrap_replicates.csv
ARD_bootstrap_CI.csv
```

## Statistical definitions

### ASIR
For age group `i`, let `d_i` be the number of events, `PY_i` the observed person-years, and `w_i` the Segi-Doll standard population weight. The directly standardized rate is:

```text
ASIR = sum(w_i * d_i / PY_i) / sum(w_i) * 100,000
```

The ASIR confidence interval is calculated with the gamma-based method from `epitools::ageadjust.direct()`.

### IRR

Age-split person-time is analyzed using a Poisson model with log person-years as the offset. Follow-up is split at attained ages 1 through 100, using the outcome-specific end age. 

```r
event ~ sex + ns(agestart, df = 4) + birth_cohort + offset(log(pyears))
```

With women as the reference level, the exponentiated male coefficient is the male-to-female IRR. ASIRs refer only to ages 40 through 84.

### ASIR difference

The absolute sex difference is defined as:

```text
ASIR difference = ASIR_male - ASIR_female
```

Its confidence interval uses a normal approximation and the sum of the independent male and female ASIR variances.

### PAF

Five-year PAFs are obtained from sex-specific Cox models using `AF::AFcoxph()`. The models contain all risk factors specified for the endpoint, together with the original adjustment variables. The PAF P value is calculated from the PAF estimate and the variance supplied by the AF object.

Joint PAFs use `AFcoxph_joint()`, which simultaneously sets all specified risk factors to zero while retaining the observed values of other covariates.

### ARD and AP

For a risk factor or the joint counterfactual:

```text
ARD = ASIR_male * PAF_male - ASIR_female * PAF_female
AP  = ARD / (ASIR_male - ASIR_female) * 100%
```

Single-factor ARDs are not additive. The joint ARD is calculated separately from the joint PAF. The AP is retained as a point estimate and is calculated as `ARD / (ASIR_male - ASIR_female)`.

## Bootstrap

The demonstration uses a small number of nonparametric bootstrap replicates to keep the synthetic demonstration practical to run. Within each replicate, males and females are resampled with replacement separately; ASIRs and Cox-based PAFs are recalculated; and ARDs are derived before taking the 2.5th and 97.5th percentiles.
