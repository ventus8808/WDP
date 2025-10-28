**Title:** Cancer Mortality Associated with County-Level Cumulative Environmental Quality

**ABSTRACT**

Background: Environmental exposures are complex and multidimensional and may contribute to cancer mortality beyond individual-level risk factors. The Environmental Quality Index (EQI) provides a composite, county-level measure of environmental burden across air, water, land, built, and sociodemographic domains.

Methods: We examined associations between county-level EQI (2000–2005 and 2006–2010) and age-adjusted mortality rates (AAMR) for all-site and site-specific cancers using four lag structures aligning exposures and outcomes (5- and 10-year lags). County-level AAMR data (ICD-10 C00–C97) were obtained from CDC WONDER (NCHS). We fitted Bayesian hierarchical mixed-effects models that accommodated interval-censored mortality data, included county-level smoking rate as a covariate, and incorporated state-level random intercepts. We compared domain-specific models and performed cluster analyses to identify environmental-regime patterns. Primary effect measures were mortality rate differences (MRD, per 100,000 population) comparing EQI quintiles (Q2–Q5 vs Q1).

Results: Across lag specifications, counties in the highest EQI quintile (Q5, poorest quality) experienced higher all-site AAMR than counties in Q1 (best quality); nationally the 5-year lag MRD was 12.12 (95% credible interval [CrI]: 9.80–14.58). Associations were strongest for the air and built domains and for cancers of the respiratory and digestive systems. Cluster analyses identified three environmental-regime clusters (high-, mixed-, low-quality); counties in the low-quality regime showed the largest MRDs (MRD = 22.51; 95% CrI: 17.62–27.68). Associations were heterogeneous by Rural–Urban Continuum Code (RUCC), climate zone, sex, and race.

Conclusions: County-level cumulative environmental degradation, as summarized by the EQI, was associated with elevated cancer mortality in the contiguous United States, with strongest signals from air and built-environment domains and in multi-domain degraded regimes. These ecological findings underscore the need for targeted environmental interventions and further individual-level and mechanistic studies to probe causal pathways.

**INTRODUCTION**

Cancer remains a leading cause of premature mortality worldwide and in the United States. Although genetic susceptibility and individual behaviors such as tobacco use are important determinants of cancer risk, involuntary environmental exposures also contribute to population-level cancer burden. Environmental exposures are multifaceted and often co-occur across media (air, water, land) and social and built environments. Composite indices such as the Environmental Quality Index (EQI) facilitate evaluation of cumulative environmental burdens at the county level. Prior studies using EQI reported associations with cancer incidence and mortality; however, questions remain about domain-specific contributions, spatial heterogeneity, and temporal alignment between exposure measurement and mortality outcomes.

This study evaluated relationships between county-level EQI and age-adjusted cancer mortality (AAMR) across multiple lag structures, domains, cancer sites, and stratified contexts. We used Bayesian hierarchical models that explicitly accounted for the interval-censored nature of mortality data from CDC WONDER and controlled for key county-level covariates.

**METHODS AND MATERIALS**

**Study design and lag structures**

We specified four exposure–outcome alignment structures: (1) EQI 2000–2005 → AAMR 2006–2010 (5-year lag), (2) EQI 2000–2005 → AAMR 2011–2015 (10-year lag), (3) EQI 2006–2010 → AAMR 2011–2015 (5-year lag), and (4) EQI 2006–2010 → AAMR 2016–2020 (10-year lag).

**Study population and outcome**

County-level AAMR for all-site and site-specific cancers (ICD-10 C00–C97) were retrieved from CDC WONDER (NCHS). Analyses included counties in the contiguous United States; Alaska and Hawaii were excluded due to spatial discontinuity. To preserve information where mortality counts were suppressed for confidentiality, we modeled AAMR as interval-censored where appropriate (exact zeros for counties with zero deaths; intervals for counties with 1–9 deaths; point estimates for counties with ≥10 deaths).

**Environmental Quality Index**

We used the U.S. EPA EQI for 2000–2005 and 2006–2010. The EQI was constructed via principal component analysis (PCA): domain-specific PCA reduced input variables to a single component per domain (air, water, land, built, sociodemographic), and a second-stage PCA combined domain indices into an overall EQI. Higher EQI values indicate poorer environmental quality. Domain and overall EQI datasets were used as published by EPA.

**Covariates**

Analyses adjusted for county-level smoking prevalence, Rural–Urban Continuum Codes (RUCC, condensed into four categories), Köppen–Geiger climate zones (tropical, dry, temperate, continental, polar; tropical and polar excluded when sample sizes were insufficient), and state-level random intercepts. RUCC and climate stratifications were used in subgroup analyses.

**Statistical models**

We fitted Bayesian hierarchical mixed-effects models to estimate associations between EQI quintiles (Q1–Q5) and AAMR while accounting for interval censoring and state-level clustering. Exact observations contributed to the likelihood via the normal density; interval-censored observations contributed via differences of cumulative distribution functions. Primary models included county-level smoking rate and state random intercepts. We reported MRDs (differences in AAMR per 100,000) and 95% credible intervals from posterior distributions. Analyses were implemented in Stan and executed from R.

We fitted complementary model frameworks: (1) overall EQI quintile models; (2) multi-domain models including quintiles for each domain simultaneously; and (3) delta models examining categorized changes in EQI (Improved, Stable, Worsened) and associated AAMR changes.

**Cluster analysis**

To identify environmental-regime patterns, we clustered counties on domain-specific EQI profiles. Algorithms compared included K-means, Gaussian Mixture Model (GMM), Birch, agglomerative hierarchical clustering, and spectral clustering. Final selection prioritized silhouette score and spatial interpretability.

**RESULTS**

**Population description**

During 2006–2010, the leading causes of cancer mortality across U.S. counties were lung, colorectal, breast, pancreatic, and prostate cancers (**Table 1**). Overall, all-site AAMR was lowest in metropolitan urbanized counties (RUCC1: 172.32 ± 0.12 per 100,000) and highest in less urbanized counties (RUCC3: 186.33 ± 0.36 per 100,000), indicating a rural–urban gradient. AAMR data for subsequent periods are provided in **Supplementary Table 1**.

**All-site cancer mortality**

Nationally, poorer cumulative environmental quality (higher EQI quintiles) was associated with higher all-site AAMR and showed dose–response patterns (**Figure 1A**, **Figure 1B**). Comparing Q5 to Q1, the 5-year lag MRD was 12.12 (95% CrI: 9.80–14.58). The air and built domains exhibited the most consistent positive associations; poorer air domain values were associated with MRD = 13.14 (95% CrI: 10.66–15.69) and the built domain with MRD = 6.77 (95% CrI: 4.33–9.14). Associations for water, land, and sociodemographic domains were less consistent and often included the null.

**Site-specific cancer mortality**

Domain-specific associations varied by cancer site (**Figure 2**). Lung cancer showed the strongest association with the air domain (Q5 vs Q1 MRD = 7.77; 95% CrI: 6.14–9.42). Colorectal cancer associated with both the air (MRD = 2.33; 95% CrI: 1.46–3.15) and land (MRD = 1.73; 95% CrI: 0.97–2.50) domains. Breast and pancreatic cancers also showed positive associations with air and built domains. Organ-system level findings are summarized in **Table 2**.

**Environmental-regime clusters**

K-means clustering (K = 3) yielded three interpretable regimes with geographic coherence (**Figure 3**; **Supplementary Figure S1**). The high-quality regime (Cluster 0) showed favorable domain values; the mixed-quality regime (Cluster 1) had heterogenous domain profiles with relatively better water quality; the low-quality regime (Cluster 2) exhibited multi-domain degradation. MRDs were largest in the low-quality regime (MRD = 22.51; 95% CrI: 17.62–27.68; **Supplementary Table S2**).

**Dynamic analyses and effect modification**

Changes in EQI were associated with corresponding changes in AAMR: environmental improvement was associated with decreases in mortality (5-year lag MRD = −2.20; 95% CrI: −3.67 to −0.33), while deterioration tended toward increases (MRD positive, marginal significance). Associations varied by regime, RUCC, climate zone, sex, and race. Sex-stratified analyses showed positive associations for both males and females, with slightly larger estimates among females. Race-stratified analyses revealed larger effect estimates among Black or African American populations compared with White populations; associations among smaller racial groups were imprecise (**Supplementary Table S3**).

**DISCUSSION**

Principal findings: In this national, county-level analysis, higher cumulative environmental burden (EQI) was associated with greater age-adjusted cancer mortality across multiple lag structures. The strongest and most consistent domain signals originated from the air and built-environment domains, and counties with multi-domain degradation experienced the largest mortality differences.

Interpretation: The positive associations with the air domain are consistent with established links between air pollution and cancer mortality, particularly for respiratory and lung cancers. The built-environment domain likely captures structural determinants related to access to services, housing quality, and neighborhood resources that shape health behaviors, exposure pathways, and healthcare access; these factors plausibly influence cancer screening, stage at diagnosis, and survival. The larger MRDs in multi-domain degraded regimes suggest cumulative or synergistic effects of co-occurring environmental stressors.

Strengths: Key strengths include use of a comprehensive, multi-domain EQI; explicit modeling of interval-censored mortality data to retain information from counties with suppressed counts; hierarchical models that accounted for state-level clustering; and thorough stratified and domain-specific analyses. The use of multiple exposure–outcome lag structures assessed temporal robustness.

Limitations: Several limitations warrant careful interpretation. First, the ecological design precludes causal inference at the individual level and is susceptible to ecological fallacy. Second, EQI is a composite index designed for comparability but does not identify specific chemical agents or mechanisms; domain contributions reflect aggregated signals. Third, residual confounding by unmeasured county-level factors (e.g., healthcare system differences, occupational exposures) may persist despite adjustment. Fourth, latency between exposure and cancer mortality is variable and often long; our lag choices are pragmatic but cannot fully capture individual-level exposure histories. Fifth, spatial dependence and regional heterogeneity may influence estimates; although state random effects mitigate some clustering, further spatial modeling may refine inference.

Sensitivity and robustness: We recommend (and performed sensitivity checks where feasible) assessing alternative EQI constructions (e.g., domain weightings), evaluating additional lag windows, testing spatial autoregressive specifications, and applying multiple-testing corrections for site-specific analyses. Where sample sizes are small (rare cancers or small racial subgroups), estimates are imprecise and should be interpreted cautiously.

Implications and policy relevance: Findings underscore population-level associations between cumulative environmental burden and cancer mortality. Interventions targeting air quality improvements and structural determinants captured in the built domain may yield public-health benefits. Identification of high-burden, low-quality regimes can inform prioritized remediation and surveillance.

Future directions: Future research should integrate individual-level cohort data, refined exposure models (e.g., spatiotemporal pollution surfaces), and mechanistic studies to evaluate causal pathways and latency. Evaluations of policy interventions (natural experiments) and more granular spatial modeling (e.g., within-county heterogeneity) are needed.

**CONCLUSION**

Counties with poorer cumulative environmental quality, as measured by the EQI, experienced higher age-adjusted cancer mortality in the contiguous United States. Air pollution and built-environment factors were the most consistent domain-level contributors, and multi-domain degraded regimes exhibited the largest mortality differences. These results highlight the relevance of cumulative environmental assessments for public-health planning and the potential value of targeted environmental and structural interventions to reduce cancer mortality burden.

**ABBREVIATIONS**

- AAMR: Age-Adjusted Mortality Rate  
- ACS: American Cancer Society  
- CDC WONDER: Centers for Disease Control and Prevention’s WONDER system  
- EQI: Environmental Quality Index  
- EPA: U.S. Environmental Protection Agency  
- GMM: Gaussian Mixture Model  
- MRD: Mortality Rate Difference  
- NCHS: National Center for Health Statistics  
- PCA: Principal Component Analysis  
- RUCC: Rural–Urban Continuum Codes

**DATA AVAILABILITY AND SOFTWARE**

EQI data are publicly available from the U.S. EPA. Mortality data were obtained from CDC WONDER (NCHS). Analyses were conducted in R (version 4.3.3) and Stan for Bayesian modeling. Code and processed outputs supporting the analyses are archived in the project repository to support reproducibility.

**ACKNOWLEDGMENTS, FUNDING, AND CONFLICTS OF INTEREST**

(Insert funding sources, acknowledgments, and conflict-of-interest disclosures here per journal requirements.)

**REFERENCES**

(References cited numerically in the text should be listed here in the journal's required format.)

**FIGURES AND TABLES (referenced in text)**

- **Figure 1A/B**: National MRD estimates by EQI quintile (5- and 10-year lags).  
- **Figure 2**: Site-specific cancer MRD estimates by domain.  
- **Figure 3**: Map of county cluster membership (environmental regimes).  
- **Figure 4**: Dynamic analysis: changes in EQI vs changes in AAMR.  
- **Table 1**: Descriptive statistics of counties and leading causes of cancer mortality.  
- **Table 2**: Organ-system level MRD estimates by lag period.  
- **Supplementary Table 1**: AAMR by period (2011–2015; 2016–2020).  
- **Supplementary Table S2**: Cluster-stratified MRD results.  
- **Supplementary Table S3**: Sex- and race-stratified MRD results.  
- **Supplementary Figures S1–S5**: Clustering diagnostics and alternate algorithm maps.
