# ANAL-AAMR Research Plan

## Study Goal

Evaluate whether county-level artificial nighttime light (ANAL) exposure is associated with later cancer and neurodegenerative disease (NDD) mortality in the United States, and whether associations differ by social vulnerability, urbanicity, and other contextual strata.

The design should use the annual ANAL panel without downloading new data. The annual panel will remain the raw exposure source, but analysis will use period-averaged exposure windows aligned to AAMR windows.

## Data Sources

- ANAL exposure: `Data/Original/ANAL/ntl_county_panel_2000_2025.csv`
  - County-year panel, 2000-2025.
  - Main exposure metric: `popw_mean_rad`.
  - Sensitivity metrics: `mean_rad`, `sol`, `lit_area_km2`, `log1p(popw_mean_rad)`.
- AAMR outcomes: existing CDC triangulation AAMR files and processed long tables.
- SVI: existing `Data/Processed/df_SVI.csv` and SVI trajectory class A/B/C/D.
- Existing contextual strata: RUCC, Census region, climate zone, economic type, homeownership tertile, Cluster_EQI, Cluster_NLCD, and available demographic strata.

## Primary Outcomes

Primary outcomes:

- Overall cancer: `C00_C97`
- Overall NDD: `G20_G30_G12.2_F01_F03`

Secondary outcomes:

- Cancer: `C34`, `C50`, `C61`, `C18_C21`, `C22`
- NDD: `G30_F01_F03`, `G20`, `G12.2`

Outcomes that are unavailable in `2021-2024` should not be emphasized in pooled models involving the most recent AAMR period.

## Exposure Definition

Primary exposure:

- Period-average `popw_mean_rad`.
- Categorize into period-specific quintiles, Q1-Q5.
- Q1 is the reference group.

Period-specific quintiles are preferred because national nighttime light levels change over time. The interpretation is:

> Within a given exposure period, counties in higher ANAL quintiles are compared with counties in the lowest ANAL quintile from the same period.

Sensitivity exposure definitions:

- Global quintiles across all exposure windows.
- Tertiles, especially for joint ANAL-SVI models.
- Continuous z-score of `log1p(popw_mean_rad)`.
- Alternative metrics: `mean_rad`, `sol`, `lit_area_km2`.

## Lag Design

Use exact calendar-shift exposure windows. Each exposure year maps to the AAMR year the specified number of years later.

| Lag | ANAL exposure window | AAMR window | Role |
|---:|---|---|---|
| 5 | 2001-2005 | 2006-2010 | Lag sensitivity |
| 5 | 2006-2010 | 2011-2015 | Lag sensitivity |
| 5 | 2011-2015 | 2016-2020 | Lag sensitivity |
| 5 | 2016-2019 | 2021-2024 | Lag sensitivity |
| 10 | 2001-2005 | 2011-2015 | Primary |
| 10 | 2006-2010 | 2016-2020 | Primary |
| 10 | 2011-2014 | 2021-2024 | Primary, with pandemic-era sensitivity |
| 15 | 2001-2005 | 2016-2020 | Lag sensitivity |
| 15 | 2006-2009 | 2021-2024 | Lag sensitivity |

The primary manuscript should emphasize the pooled 10-year lag because it is concise, biologically plausible for cancer/NDD mortality, and gives one clear conclusion.

## Main Analysis

Run three pooled lag-specific models:

- Pooled 5-year lag model.
- Pooled 10-year lag model.
- Pooled 15-year lag model.

The main result should be the pooled 10-year lag model.

Recommended pooled model:

```text
AAMR interval ~ ANAL_quintile + matched_pair + disease-specific covariates
              + state random intercept + county random intercept
```

If county random intercept is too expensive or unstable:

```text
AAMR interval ~ ANAL_quintile + matched_pair + disease-specific covariates
              + state random intercept
```

Also run single-pair models for each exposure-outcome window. These are not the main conclusion; they show temporal consistency.

## Covariate Strategy

Use disease-specific covariate sets consistent with the current SVI pipeline:

- Cancer: `Smoking_rate + Uninsured_rate`
- NDD: `Physical_Activities_rate + Uninsured_rate`

Candidate sensitivity covariates:

- `Obesity_rate`
- `Physician_Density_per100k`
- `Diabetes_Prevalence_rate`
- `Forest_Coverage`
- `RUCC`
- `SVI`

Avoid adjusting for both highly related urbanicity measures and stratifying by urbanicity as the only interpretation. For example, if RUCC is a stratification variable, also present a non-RUCC-adjusted main model.

## SVI Combination

SVI should be used in two ways.

First, as an adjustment variable:

```text
AAMR interval ~ ANAL_quintile + SVI + matched_pair + covariates
              + random effects
```

Second, as an effect modifier:

```text
AAMR interval ~ ANAL_tertile * SVI_class + matched_pair + covariates
              + random effects
```

Use `ANAL_tertile x SVI A/B/C/D` as the primary joint exposure. This creates 12 groups and is more readable than `ANAL_quintile x SVI`, which creates 20 groups.

Reference group:

- Low ANAL tertile + SVI class A.

## Stratified Analysis

Stratification should be targeted, not exhaustive in the main manuscript.

Primary stratification:

1. RUCC, because ANAL is strongly related to urbanicity.
2. SVI class, preferably through interaction/joint exposure.
3. Census region, to check geographic consistency.

Secondary or supplement stratification:

- Climate zone.
- Economic type.
- Homeownership tertile.
- Cluster_EQI.
- Cluster_NLCD.
- Available sex/race demographic strata from `df_SVI_Stratified.csv`.

For stratified analysis, start with:

- Overall cancer and overall NDD only.
- Pooled 10-year lag only.
- `popw_mean_rad` quintile or tertile.

Expand to secondary outcomes only after the primary stratified results are stable and interpretable.

## Test 1: Covariate Combination Test

Purpose:

- Determine whether the ANAL-AAMR association is robust to different confounder sets.
- Avoid overclaiming from a single adjustment model.

Primary setting:

- Use pooled 10-year lag.
- Use overall cancer and overall NDD first.
- Exposure: period-specific `popw_mean_rad` quintile.

Model sequence:

| Model | Adjustment set | Purpose |
|---|---|---|
| M0 | ANAL + matched_pair + random effects | Crude temporal/context model |
| M1 | M0 + disease-specific core covariates | Primary model |
| M2 | M1 + SVI | Social vulnerability adjustment |
| M3 | M1 + RUCC | Urbanicity adjustment |
| M4 | M1 + SVI + RUCC | Combined social/urban adjustment |
| M5 | M1 + OB + PD + DB + FC | Expanded health/environmental covariates |
| M6 | M1 + SVI + RUCC + OB + PD + DB + FC | Fully adjusted sensitivity |

Abbreviations:

- OB: `Obesity_rate`
- PD: `Physician_Density_per100k`
- DB: `Diabetes_Prevalence_rate`
- FC: `Forest_Coverage`

Evaluation:

- Compare Q5 vs Q1 direction, magnitude, credible interval, and posterior tail probability.
- Check monotonicity across Q2-Q5.
- Check R-hat and ESS.
- Identify whether the association is stable, attenuated by SVI/RUCC, or only present in crude models.

Optional extended test:

- All subsets of candidate covariates, similar to the existing SVI sensitivity-combination scripts.
- Use only for overall cancer and overall NDD because the full combination grid can become large.

## Test 2: Lag Test

Purpose:

- Determine whether the association is most consistent at 5, 10, or 15 years.
- Avoid presenting all nine windows as equally primary.

Primary lag test:

- Run pooled lag-specific models:
  - pooled 5-year lag
  - pooled 10-year lag
  - pooled 15-year lag

Secondary lag test:

- Run all single-pair models:
  - 4 single 5-year pairs
  - 3 single 10-year pairs
  - 2 single 15-year pairs

Evaluation:

- Main evidence: pooled 10-year lag.
- Robustness: similar direction across 5-, 10-, and 15-year pooled models.
- Consistency: single-pair estimates should not be driven by one period only.
- Recent-period check: repeat lag models excluding pairs ending in `2021-2024`.

Recommended interpretation:

- If 10-year pooled results are strongest and consistent, use 10-year lag as the primary manuscript analysis.
- If all lags show similar direction, report broad lag robustness.
- If only 5-year lag is present, interpret as short-term mortality/progression association rather than long-latency etiology.
- If only 15-year lag is present, interpret as longer-latency association, especially for cancer.

## Sensitivity Analyses

Exposure sensitivity:

- Replace `popw_mean_rad` with `mean_rad`, `sol`, and `lit_area_km2`.
- Use continuous `log1p(popw_mean_rad)` z-score.
- Compare period-specific vs global quintiles.

Time-period sensitivity:

- Exclude all pairs ending in `2021-2024`.
- Exclude 2024 if concerned about the 2023 to 2024 source transition.
- Compare pooled 10-year model with and without `2011-2014 -> 2021-2024`.

Model sensitivity:

- State random intercept only.
- State + county random intercept.
- Add matched-pair fixed effects in all pooled models.

Outcome sensitivity:

- Repeat primary pooled 10-year model for secondary cancer and NDD outcomes.
- Avoid strong conclusions for sparse or incomplete outcome-period combinations.

## Suggested Reporting Structure

Main manuscript:

1. Pooled 10-year ANAL-AAMR result for overall cancer and overall NDD.
2. Pooled 5-, 10-, and 15-year comparison figure.
3. ANAL-SVI joint exposure result for overall cancer and NDD.
4. RUCC-stratified pooled 10-year result.

Supplement:

- All nine single-pair lag models.
- Covariate combination test.
- Alternative exposure metrics.
- Global vs period-specific quintiles.
- Secondary outcomes.
- Additional strata.

## Recommended Workflow

1. Build ANAL-AAMR analysis table with the nine lag-pair windows.
2. Run pooled 10-year main models for overall cancer and NDD.
3. Run single-pair 10-year models to check consistency.
4. Run pooled 5-year and 15-year lag tests.
5. Run covariate combination test on pooled 10-year models.
6. Run ANAL-SVI joint model.
7. Run RUCC-stratified pooled 10-year model.
8. Expand to secondary outcomes and supplementary strata.

