# Publication-Quality Ridgeline Figures

**Generated:** 2024-12-11  
**Design Standard:** Nature / The Lancet  
**Resolution:** 600 DPI (publication-ready)

---

## 📊 Figure Inventory

### Main Figures (Manuscript Body)

#### 1. `Fig_Domains_ByLag_Publication.png` (1.1 MB)
- **Layout:** 6 domains (rows) × 3 lag periods (columns)
- **Content:** Comprehensive dose-response across all EQI domains
- **Domains:** Overall EQI, Air, Water, Land, Built, Social
- **Lag Periods:** 5, 10, 15 years
- **Color Palette:** Warm gradient (light orange → deep red)
- **Recommended Use:** Main text Figure 2 or 3 (primary results)
- **Key Features:**
  - Shows domain-specific effects over time
  - Identifies which environmental factors drive mortality
  - Reveals temporal heterogeneity in effects

#### 2. `Fig_Overall_Evolution_Publication.png` (588 KB)
- **Layout:** Q2-Q5 quintiles × 3 lag periods (stacked)
- **Content:** Temporal dynamics of aggregate EQI effects
- **Color Palette:** Cool gradient (light blue → navy)
- **Recommended Use:** Main text Figure 1 or Supplementary
- **Key Features:**
  - Simplified view of overall environmental burden
  - Clear dose-response gradient (Q2→Q5)
  - Temporal evolution of effect sizes

### Supplementary Figures (Detailed Views)

#### 3. `Fig_Lag5_Domains_Publication.png` (452 KB)
- **Period:** Short-term effects (5-year lag, 2006-2010 AAMR)
- **Layout:** 6 domains stacked vertically
- **Use:** Supplementary material for detailed inspection

#### 4. `Fig_Lag10_Domains_Publication.png` (451 KB)
- **Period:** Medium-term effects (10-year lag, 2011-2015 AAMR)
- **Layout:** 6 domains stacked vertically
- **Use:** Supplementary material for detailed inspection

#### 5. `Fig_Lag15_Domains_Publication.png` (449 KB)
- **Period:** Long-term effects (15-year lag, 2016-2020 AAMR)
- **Layout:** 6 domains stacked vertically
- **Use:** Supplementary material for detailed inspection

---

## 🎨 Design Philosophy

### Aesthetic Minimalism
- **Single-color gradients:** Q2 → Q5 represented by intensity (light to dark)
- **Clean typography:** Helvetica, Nature/Lancet font sizes
- **Minimal grid lines:** Only major x-axis gridlines retained
- **White background:** Maximum contrast and clarity
- **No decorative elements:** Data-ink ratio maximized

### Information Maximalism
- **Bayesian probabilities:** P(effect > 0) annotated on each distribution
- **Zero reference line:** Dashed vertical line at MRD = 0 (null effect)
- **Posterior quantiles:** White lines show median (and 2.5%/97.5% in some views)
- **Dose-response gradient:** Visual metaphor of "risk waves" moving rightward
- **Complete uncertainty:** Full posterior distributions (not just point estimates)

---

## 📈 Key Visual Elements

### Color Gradients

**Warm Palette** (Domain-specific figures):
- Q2: `#FED9B7` (light peach)
- Q3: `#F4A261` (soft orange)
- Q4: `#E76F51` (coral)
- Q5: `#D62828` (crimson)
- Q5+: `#9B2226` (deep red)

**Cool Palette** (Overall evolution):
- Q2: `#BDE0FE` (pale blue)
- Q3: `#A2D2FF` (sky blue)
- Q4: `#6A9FD4` (cornflower)
- Q5: `#4A7BA7` (steel blue)
- Q5+: `#2C5F8D` (navy)

### Annotations

**Posterior Probabilities:**
- `P>0.999` — virtually certain positive effect
- `P=0.XX` — standard probability format
- `P<0.001` — virtually certain negative/null effect
- **Interpretation:** Probability that mortality risk is elevated (vs. Q1 reference)

**Reference Lines:**
- Dashed black line at x = 0 (MRD = 0) indicates null effect
- Area to right (positive): harmful effect
- Area to left (negative): protective effect

---

## 📊 Data Summary

### Model Specifications
- **Method:** Bayesian interval-censored spatiotemporal mixed-effects model
- **Outcome:** Age-Adjusted Mortality Rate (AAMR) for all cancers (C00-C97)
- **Exposure:** Environmental Quality Index (EQI) quintiles
- **Reference:** Q1 (lowest EQI, best environmental quality)
- **Spatial Units:** 3,106 US counties
- **Temporal Coverage:** EQI 2000-2005 → AAMR 2006-2020 (lagged)

### Key Findings (from figures)

**Air Domain (strongest effects):**
- Lag 5: Q5 effect = +13.5 [+10.6, +16.3], P > 0.999
- Lag 10: Q5 effect = +12.3 [+9.4, +15.2], P > 0.999
- Consistent strong dose-response across all lags

**Water Domain (weak/uncertain):**
- Effects generally small and uncertain
- Most quintiles cross zero (null effect line)
- P(effect > 0) ≈ 0.6-0.9 (moderate probability)

**Social Domain (protective effects):**
- Negative effects (protective) across multiple quintiles
- Lag 5: Q3 effect = -5.8 [-7.9, -3.5], P < 0.001
- Unexpected inverse association

**Built Domain (moderate effects):**
- Consistent positive effects across lags
- Q5 effects range from +6 to +8 AAMR units

---

## 🔬 Interpretation Guide

### Reading Ridgeline Plots

1. **Y-axis:** EQI quintiles (Q2-Q5 vs. Q1 reference) or domains
2. **X-axis:** Effect on Age-Adjusted Mortality Rate (deaths per 100,000)
3. **Distribution shape:** Full posterior uncertainty from Bayesian model
4. **White vertical lines:** Posterior median (50th percentile)
5. **Shading intensity:** Darker = higher exposure level (worse EQI)

### Statistical Interpretation

- **Positive effects:** Environmental degradation → increased mortality
- **Negative effects:** Better environmental indicators → increased mortality (paradoxical)
- **Wide distributions:** High uncertainty, need more data
- **Narrow distributions:** Precise estimates, strong evidence
- **P-values:** Bayesian posterior probabilities (NOT frequentist p-values)

### Clinical/Public Health Interpretation

**Example (Air, Q5, Lag 5):**
> "Counties in the worst air quality quintile (Q5) experienced 13.5 additional 
> cancer deaths per 100,000 population (95% CrI: 10.6-16.3) compared to the best 
> quintile (Q1), with posterior probability > 0.999 that the effect is positive."

---

## 📝 Recommended Figure Captions

### Main Figure 1 (Domains × Lag)
```
Environmental Quality Index domain-specific effects on cancer mortality across 
lag periods. Posterior distributions show quintile effects (Q2-Q5 vs. reference Q1) 
from Bayesian interval-censored spatiotemporal models. Gradient shading represents 
exposure intensity (lighter = lower EQI quintile). White lines indicate posterior 
median. Dashed vertical line marks null effect (MRD = 0). Italic annotations show 
posterior probability of positive effect. Air domain shows strongest and most 
consistent associations, while Social domain exhibits unexpected protective effects. 
Data: EQI 2000-2005, AAMR 2006-2020, N=3,106 US counties.
```

### Main Figure 2 (Overall Evolution)
```
Temporal evolution of overall Environmental Quality Index effects on cancer 
mortality. Dose-response gradient (Q2→Q5, light to dark blue) shows consistent 
positive associations across all lag periods, with effect sizes declining over 
time (short-term: Lag 5 > medium-term: Lag 10 > long-term: Lag 15). White lines 
show posterior median and 95% credible intervals. Bayesian posterior probabilities 
(P) indicate near-certainty of positive effects across most quintiles.
```

---

## 🛠️ Technical Details

### Generation
- **Script:** `Code/Visualization/RidgeLine_Publication.R`
- **Dependencies:** R packages: `ggplot2`, `ggridges`, `dplyr`, `tidyr`, `scales`
- **Input:** `.rds` files from `Code/brms/cmdstan_main_ridgeline.R`
- **Runtime:** ~2-3 minutes per figure (600 DPI rendering)

### Specifications
- **Resolution:** 600 DPI (Nature/Lancet standard)
- **Format:** PNG with white background
- **Color Space:** sRGB
- **Font:** Helvetica (or system sans-serif)
- **Dimensions:** 
  - Main (Domains×Lag): 12" × 14"
  - Evolution: 10" × 10"
  - Supplementary: 10" × 8"

### Reproducibility
To regenerate all figures:
```bash
bash Code/Visualization/generate_publication_figures.sh
```

Or individual figures:
```bash
# Main figure: all domains by lag
Rscript Code/Visualization/RidgeLine_Publication.R \
  --mode domains_by_lag --palette warm

# Overall temporal evolution
Rscript Code/Visualization/RidgeLine_Publication.R \
  --mode overall_evolution --palette cool

# Supplementary: single lag detail
Rscript Code/Visualization/RidgeLine_Publication.R \
  --mode domains_single --lag 10 --palette warm
```

---

## 📋 Checklist for Manuscript Submission

- [x] All figures generated at 600 DPI
- [x] Posterior probabilities annotated
- [x] Zero reference lines present
- [x] Color gradients represent dose-response
- [x] Minimalist Nature/Lancet aesthetics applied
- [x] Font sizes appropriate for print (8-12 pt)
- [ ] Figures reviewed by co-authors
- [ ] Color-blind friendly palette confirmed (warm = distinguishable in grayscale)
- [ ] Figure captions drafted and reviewed
- [ ] Supplementary figure legends prepared
- [ ] High-resolution files backed up

---

## 📚 Related Files

- **Data Processing:** `Code/brms/cmdstan_main_ridgeline.R`
- **Visualization Code:** `Code/Visualization/RidgeLine_Publication.R`
- **Batch Generator:** `Code/Visualization/generate_publication_figures.sh`
- **Original Test Script:** `Code/Visualization/RidgeLine_Test.R`
- **Model Outputs:** `Result/Ridgeline/*.rds`

---

## 🔄 Version History

- **v1.0** (2024-12-11): Initial publication-quality figures
  - Implemented monochromatic gradients
  - Added Bayesian posterior probabilities
  - Nature/Lancet minimalist theme
  - 600 DPI resolution

---

## 📧 Contact

For questions about figure generation or interpretation:
- See: `.github/copilot-instructions.md`
- Repository: WDP (Environmental Quality and Cancer Mortality)

---

**End of Documentation**