# County Economic Typology × EQI Analysis Report

## Summary Statistics

### Mean EQI Scores by County Economic Type

| Economic Type | N | Total EQI | Air | Water | Land | Built | Social |
|---------------|---|-----------|-----|-------|------|-------|--------|
| Farming | 440 | -0.891 | -1.300 | -0.245 | 0.180 | -0.848 | -0.349 |
| Mining | 128 | -0.574 | -0.555 | 0.212 | -0.261 | -0.187 | -0.731 |
| Manufacturing | 905 | 0.210 | 0.423 | 0.044 | 0.030 | -0.007 | 0.126 |
| Government | 381 | -0.091 | 0.082 | -0.073 | -0.280 | 0.136 | -0.202 |
| Services | 339 | 0.811 | 0.611 | 0.074 | 0.056 | 0.827 | 0.700 |
| Nonspecialized | 948 | 0.037 | 0.024 | 0.046 | 0.015 | 0.075 | -0.029 |

## Statistical Tests

### ANOVA Results

| Domain | F-statistic | P-value | η² | Significant |
|--------|-------------|---------|-----|-------------|
| EQI | 166.68 | 0.000000 | 0.210 | Yes |
| EQI_Air | 326.20 | 0.000000 | 0.342 | Yes |
| EQI_Water | 8.04 | 0.000000 | 0.013 | Yes |
| EQI_Land | 11.16 | 0.000000 | 0.017 | Yes |
| EQI_Built | 137.65 | 0.000000 | 0.180 | Yes |
| EQI_Social | 70.74 | 0.000000 | 0.101 | Yes |

### Key Findings

- **EQI**: Significant differences across county types (F=166.68, p<0.001, large effect size η²=0.210)
- **EQI_Air**: Significant differences across county types (F=326.20, p<0.001, large effect size η²=0.342)
- **EQI_Water**: Significant differences across county types (F=8.04, p<0.001, small effect size η²=0.013)
- **EQI_Land**: Significant differences across county types (F=11.16, p<0.001, small effect size η²=0.017)
- **EQI_Built**: Significant differences across county types (F=137.65, p<0.001, large effect size η²=0.180)
- **EQI_Social**: Significant differences across county types (F=70.74, p<0.001, medium effect size η²=0.101)

### Pairwise Comparisons (Selected Significant Results)

| Domain | Comparison | Mean Diff | Cohen's d | P-value |
|--------|------------|-----------|-----------|---------|
| EQI_Built | Services vs Nonspecialized | 0.751 | 0.913 | 0.000000 |
| EQI | Services vs Nonspecialized | 0.774 | 0.864 | 0.000000 |
| EQI_Social | Services vs Nonspecialized | 0.729 | 0.718 | 0.000000 |
| EQI_Air | Services vs Nonspecialized | 0.587 | 0.694 | 0.000000 |
| EQI_Air | Manufacturing vs Nonspecialized | 0.399 | 0.548 | 0.000000 |
| EQI_Social | Mining vs Government | -0.529 | -0.518 | 0.000001 |
| EQI_Air | Government vs Services | -0.528 | -0.550 | 0.000000 |
| EQI_Social | Farming vs Manufacturing | -0.475 | -0.558 | 0.000000 |
| EQI_Built | Farming vs Mining | -0.662 | -0.562 | 0.000000 |
| EQI_Social | Manufacturing vs Services | -0.573 | -0.630 | 0.000000 |
| EQI_Air | Mining vs Government | -0.638 | -0.680 | 0.000000 |
| EQI | Mining vs Nonspecialized | -0.612 | -0.688 | 0.000000 |
| EQI | Manufacturing vs Services | -0.602 | -0.690 | 0.000000 |
| EQI_Air | Mining vs Nonspecialized | -0.579 | -0.715 | 0.000000 |
| EQI_Social | Mining vs Nonspecialized | -0.702 | -0.715 | 0.000000 |
| EQI_Built | Government vs Services | -0.691 | -0.780 | 0.000000 |
| EQI_Social | Government vs Services | -0.901 | -0.845 | 0.000000 |
| EQI_Air | Farming vs Mining | -0.745 | -0.863 | 0.000000 |
| EQI | Farming vs Government | -0.801 | -0.874 | 0.000000 |
| EQI_Built | Farming vs Government | -0.984 | -0.881 | 0.000000 |

## Interpretation

This analysis examines the relationship between county economic dependency types 
(as classified by USDA ERS 2004) and environmental quality across multiple domains.

**Key Observations:**

- Counties with the **lowest** overall environmental quality: **Farming** (mean EQI = -0.891)
- Counties with the **highest** overall environmental quality: **Services** (mean EQI = 0.811)
