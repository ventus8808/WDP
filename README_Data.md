# WDP Data Sources and Processing Guide

This document records the full pipeline from Original Data → Cleaning/Processing → Analysis-ready datasets in the WDP project. As we migrate/standardize scripts, please continue to update this file with each data source’s inputs, outputs, variable definitions, and assumptions to ensure reproducibility and maintainability.

## Contents

- CDC (Outcomes and Geographic Classification)
- Pesticide (USGS PNSP)
- Socioeconomic (SAIPE / LAUS / BEA / USDA-ERS / SEER / County Adjacency)
- Environmental (CACES LUR / GEE / NLDAS)
- PCA-Derived Covariates (SVI and Climate Factors)

---

## CDC (Outcomes and Geographic Classification)

### 1) Location (Static County-Level Geography)

- **Cleaning script**: `Code/Clean/CDC_Location_Urbanization.py`
- **Input directory**: `Data/Original/ CDC WONDER/Location and Urbanization`
  - **Files used**:
    - `Location_HHS_State.csv` (HHS regions)
    - `Location_Region_Division_State.csv` (Census regions/divisions)
- **Output file**: `Data/Processed/CDC/Location.csv`
- **Granularity**: County-level (one row per county; static across years)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `County`: County name with state abbreviation
  - `HHS_Region`: HHS region (e.g., “HHS Region #1 …”)
  - `Census_Region`: Census region (e.g., Northeast/West)
  - `Census_Division`: Census division (e.g., New England)

**Use**: Static join table for any county-level panel; join via `COUNTY_FIPS`.

---

### 2) Urbanization (County × Year Panel)

- **Cleaning script**: `Code/Clean/CDC_Location_Urbanization.py`
- **Input directory**: `Data/Original/ CDC WONDER/Location and Urbanization`
  - **Files used**: `Location_County_Urbanization*.csv` (spanning 1999–2020)
- **Output file**: `Data/Processed/CDC/Urbanization.csv`
- **Granularity**: County × Year (panel)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (1999–2020)
  - `County`: County name with state abbreviation
  - `Urbanization_Code`: Code from “2013 Urbanization Code”
  - `Urbanization_Type`: Category label (e.g., Large Central Metro, NonCore)

**Processing rules**:
- For duplicate county–year records, keep the first occurrence.
- Year is coerced to integer; records with missing year are dropped.

**Use**: Yearly urbanization classification for county-level analyses; join via `COUNTY_FIPS` and `Year`.

---

### 3) Cancer Mortality (ICD Group Merge)

- **Cleaning script**: `Code/Clean/CDC_Death_Merge.py`
- **Input base directory**: `Data/Original/ CDC WONDER/`
  - **Example**: `Data/Original/ CDC WONDER/C00-C97`
- **Output file**: `Data/Processed/CDC/<ICD>.csv` (e.g., `C00-C97.csv`)
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int, filtered to 1999–2019)
  - `County`: County name with state abbreviation
  - `Sex`: Sex code
  - `Race`: Race category
  - `Age`: Ten-year age group
  - `Deaths`: Count of deaths (int)
  - `Population`: Population (int)
  - `SD`: Crude rate standard error

**Processing rules and notes**:
- Paths are resolved from `config.yaml` (`data_sources.cdc_wonder.integrity_base_dir`).
- Set `MANUAL_ICD_GROUP` at the top of the script to select the disease group directory.
- Merges all CSVs in the specified ICD folder.

---

### 4) AAMR (Age-Adjusted Mortality Rate, County × Year, Aggregated)

- **Cleaning script**: `Code/Clean/CDC_AAMR_Merge.py`
- **Input directory**: `Data/Original/CDC WONDER AAMR/{ICD}`（例如 `C81-C96`）
- **Output file**: `Data/Processed/CDC/AAMR_{ICD}.csv`（例如 `AAMR_C81-C96.csv`）
- **Granularity**: County × Year（1999–2020）
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `County`: County name with state abbreviation
  - `Deaths`: Death count (Int64，可空)
  - `Population`: Population (Int64，可空)
  - `CMR`: Crude Mortality Rate
  - `CMR_SE`: Standard error of CMR
  - `AAMR`: Age-adjusted mortality rate
  - `AAMR_SE`: Standard error of AAMR

**Processing rules**:
- Drop non-county totals and invalid FIPS（移除诸如全国合计 `00nan` 等记录）。
- Standardize FIPS to 5-digit string；`Deaths`/`Population` 保存为可空整型（缺失保留为 NaN）。
- Year coerced to integer; filter to 1999–2020.
- Deduplicate county-year by preferring non-missing SE; then non-missing deaths.

**Use**: Aggregated outcome for spatio-temporal Bayesian smoothing（策略三），在 `Code/Test/BSTM_Run.py` 中可通过 `--aggregated` 使用。

---

### 5) CDC Data Integrity Checker

- **Utility script**: `Code/Clean/CDC_Data_Integrity_Checker.py`
- **Input directory**: Configurable via `MANUAL_ICD_GROUPS` (e.g., "C00-C14")
- **Purpose**: Validates data completeness and structure across CDC WONDER files
- **Function**: Checks for required columns, data ranges, and file consistency

**Key features**:
- Validates presence of key columns (Year, County, County Code, Sex, Race, Age)
- Checks statistical columns (Deaths, Population, Crude Rate Standard Error)
- Reports year ranges and missing data patterns
- Compact summary mode for efficient reporting
- Configurable to check single or multiple ICD disease groups

**Use**: Run before `CDC_Death_Merge.py` to identify data quality issues

---

## Pesticide (USGS PNSP)

### 1) PNSP Compound Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/PNSP_Merge.py`
- **Input directory**: `Data/Original/USGS PNSP`
  - **Files used**: Individual year files (1999-2012, 2018-2019) + combined file (2013-2017)
  - **Mapping file**: `Data/Processed/Pesticide/mapping.csv`
- **Output file**: `Data/Processed/Pesticide/PNSP.csv`
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int, 1999–2019)
  - `cat{id}_min/avg/max`: Category-level estimates in kg/county/year (35 categories)
  - `chem{id}_min/avg/max`: Compound-specific estimates in kg/county/year (509 compounds)

**Processing rules and notes**:
- Paths are read from `config.yaml`.
- Merges all PNSP files from 1999-2019.
- Uses `mapping.csv` to map compound names to IDs and categories.

### 2) PNSP Pesticide Density Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/PNSP_Density.py`
- **Input files**:
  - `Data/Processed/Pesticide/PNSP.csv` (Pesticide weights)
  - `Data/Processed/Environmental/NLCD_JRC.csv` (Agricultural land area)
- **Output file**: `Data/Processed/Pesticide/PNSP_Density.csv`
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int, 1999–2019)
  - `cat{id}_min/avg/max`: Category-level estimates in kg/km²/year
  - `chem{id}_min/avg/max`: Compound-specific estimates in kg/km²/year

**Processing rules and notes**:
- Calculates pesticide density by dividing the application weight from `PNSP.csv` by the agricultural land area (`nlcd_agriculture_km2`) from `NLCD_JRC.csv`.
- This provides a measure of exposure intensity that accounts for the size of a county's agricultural footprint.
- If agricultural area is zero or missing, the resulting density is zero.

### 3) Compound Name Mapping and Utility Scripts

- **Utility script**: `Code/Clean/PNSP_Unique_Name.py`
- **Input directory**: `Data/Original/USGS PNSP` (all .txt and .csv files)
- **Output file**: `Data/Processed/Pesticide/Unique_Name.txt`
- **Purpose**: Extracts unique pesticide compound names from all PNSP files.

**File**: `Data/Processed/Pesticide/mapping.csv`
- **Contains**: 509 compounds with classification into 35 categories.
- **Used by**: `PNSP_Merge.py` for data reshaping and categorization.

---

## Socioeconomic

### 1) SAIPE Poverty and Income Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/SE_SAIPE_Merge.py`
- **Input directory**: `Data/Original/SAIPE`
- **Output file**: `Data/Processed/Socioeconomic/Poverty_Income.csv`
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `Poverty_Percent_All_Ages`: Poverty rate for all ages (%)
  - `Median_Household_Income`: Median household income (dollars)

### 2) LAUS Labor Force and Unemployment Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/SE_LAUS_Merge.py`
- **Input directory**: `Data/Original/LAUS`
- **Output file**: `Data/Processed/Socioeconomic/Unemployment.csv`
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `Labor_Force`: Total labor force (integer)
  - `Employed`: Number of employed persons (integer)
  - `Unemployed`: Number of unemployed persons (integer)
  - `Unemployment_Rate`: Unemployment rate (%)

### 3) BEA Economic Indicators Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/SE_BEA_Merge.py`
- **Input directory**: `Data/Original/BEA`
- **Output file**: `Data/Processed/Socioeconomic/GDP.csv`
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `Population`: Population count (integer)
  - `Total_GDP_10K_USD`: Total GDP in 10,000 USD (float)
  - `Per_Capita_Income`: Per capita income (integer)

### 4) USDA ERS Education Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/SE_USDA_ERS_Education_Merge.py`
- **Input directory**: `Data/Original/USDA ERS`
- **Output file**: `Data/Processed/Socioeconomic/Education.csv`
- **Granularity**: County × Year (panel; 1999–2020, interpolated)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `Less_Than_High_School_Percent`: Less than high school education (%)
  - `High_School_Only_Percent`: High school diploma only (%)
  - `Some_College_Percent`: Some college or associate degree (%)
  - `College_Plus_Percent`: Four-year college or higher (%)
  - `Rural_Urban_Continuum_Code`: Rural-urban continuum classification code
  - `Urban_Influence_Code`: Urban influence classification code

### 5) SEER Population Structure Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/SE_SEER_Population.py`
- **Input directory**: `Data/Original/SEER Population`
  - **File used**: `us.1990_2023.singleages.through89.90plus.adjusted.txt`
- **Output file**: `Data/Processed/Socioeconomic/Population_Structure.csv`
- **Granularity**: County × Year (panel; 1999–2020)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int, 1999–2020)
  - `Race`: Numeric code (1: White, 2: Black, 3: American Indian/Alaska Native, 4: Asian or Pacific Islander)
  - `Origin`: Numeric code (0: Non-Hispanic, 1: Hispanic)
  - `Sex`: Numeric code (1: Male, 2: Female)
  - `Age`: Numeric code (0-89 for single years, 90 for 90+)
  - `Population`: Population count (integer)

### 6) County Adjacency Data (Spatial Relationships)

- **Cleaning script**: `Code/Clean/County_Adjacency.py`
- **Input directory**: `Data/Original/County Shapeline`
  - **File used**: `tl_2015_us_county.shp` (county shapefile)
- **Output files**: 
  - `Data/Processed/Socioeconomic/County_Adjacency_Matrix.csv`
  - `Data/Processed/Socioeconomic/County_Adjacency_List.csv`
- **Granularity**: County spatial relationships (static)
- **Variables**:
  - **Matrix format** (`County_Adjacency_Matrix.csv`):
    - Rows and columns: County GEOIDs (5-digit FIPS codes)
    - Values: Boolean (True = counties are adjacent, False = not adjacent)
  - **Edge list format** (`County_Adjacency_List.csv`):
    - `county_from`: Origin county GEOID (string)
    - `county_to`: Destination county GEOID (string) 
    - `adjacency_weight`: Boolean adjacency indicator (True for adjacent counties)

**Processing rules and notes**:
- Uses GeoPandas to determine spatial adjacency based on touching boundaries
- Matrix is symmetric (if county A is adjacent to county B, then B is adjacent to A)
- Self-adjacency is set to False (counties are not adjacent to themselves)
- Edge list contains only True adjacency relationships (18,962 relationships total)
- Used for spatial modeling in Bayesian analysis (ICAR/BYM models)

---

## Environmental

### 1) CACES LUR Air Pollution Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/ENV_LUR_Merge.py`
- **Input directory**: `Data/Original/CACES LUR`
  - **Files used**: 
    - `1999-2019 O3.csv`
    - `1999-2020 CO SO2 NO2 PM10 PM25.csv`
- **Output file**: `Data/Processed/Environmental/Air_Pollution.csv`
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `O3`: Ozone concentration (pred_wght)
  - `CO`: Carbon Monoxide concentration (pred_wght)
  - `SO2`: Sulfur Dioxide concentration (pred_wght)
  - `NO2`: Nitrogen Dioxide concentration (pred_wght)
  - `PM10`: Particulate Matter < 10µm concentration (pred_wght)
  - `PM25`: Particulate Matter < 2.5µm concentration (pred_wght)

**Processing rules**:
- Removes geographical coordinates (lat, lon) and state abbreviations
- Converts from long format (pollutant column) to wide format (one column per pollutant)
- Filters years to 1999-2019 for consistency with other datasets

---

### 2) GEE Land Cover and Surface Water Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/ENV_GEE_Merge.py`
- **Input directory**: `Data/Original/GEE`
  - **Files used**: 
    - `county_base*.csv` (county area data)
    - `jrc_water_*.csv` (JRC surface water by year)
    - `nlcd_landuse_*.csv` (NLCD land cover by year)
- **Output file**: `Data/Processed/Environmental/NLCD_JRC.csv`
- **Granularity**: County × Year (panel; 1999–2020, interpolated)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `total_area_km2`: Total county area in square kilometers
  - `jrc_permanent_water_km2`: Area of permanent surface water (JRC)
  - `jrc_seasonal_water_km2`: Area of seasonal surface water (JRC)
  - `nlcd_forest_km2`: Forested area (NLCD)
  - `nlcd_water_km2`: Water area (NLCD)
  - `nlcd_urban_km2`: Urban/developed area (NLCD)
  - `nlcd_agriculture_km2`: Agricultural area (NLCD)
  - `nlcd_cropland_km2`: Cropland area (NLCD)
  - `nlcd_pasture_km2`: Pasture/hay area (NLCD)
  - `nlcd_wetland_km2`: Wetland area (NLCD)
  - `nlcd_shrub_km2`: Shrub/scrub area (NLCD)
  - `nlcd_grassland_km2`: Grassland/herbaceous area (NLCD)
  - `nlcd_barren_km2`: Barren land area (NLCD)

**Processing rules**:
- Merges county base data with JRC water and NLCD land cover data
- Creates complete county-year panel for 1999-2020
- Applies linear interpolation to fill missing yearly data
- Values are rounded to 4 decimal places and clipped to non-negative

**GEE Data Collection Scripts** (Google Earth Engine):
- `Code/Clean/ENV_GEE_County.py`: Exports county base data with total area
- `Code/Clean/ENV_GEE_JRC.py`: Exports JRC water data by year (1999-2020)
- `Code/Clean/ENV_GEE_NLCD.py`: Exports NLCD land cover data by available years
- `Code/Clean/ENV_GEE_utils.py`: Utility functions for GEE processing

---

### 3) NLDAS Meteorology Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/ENV_NLDAS_Meterology.py`
- **Input directories**:
  - `Data/Original/NLDAS` (NLDAS NetCDF files)
  - `Data/Original/County Shapeline/` (county boundaries)
- **Output file**: `Data/Processed/Environmental/NLDAS_{START_YEAR}_{END_YEAR}.csv`
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `year`: Calendar year (int)
  - `tas_mean_annual/DJF/MAM/JJA/SON`: Mean temperature (°C)
  - `wind_mean_annual/DJF/MAM/JJA/SON`: Mean wind speed (m/s)
  - `prcp_sum_annual/DJF/MAM/JJA/SON`: Total precipitation (mm/month)
  - `rh_mean_annual/DJF/MAM/JJA/SON`: Mean relative humidity (%)
  - `swrad_mean_annual/DJF/MAM/JJA/SON`: Mean shortwave radiation (W/m²)
  - `lwrad_mean_annual/DJF/MAM/JJA/SON`: Mean longwave radiation (W/m²)
  - `psurf_mean_annual/DJF/MAM/JJA/SON`: Mean surface pressure (kPa)
  - `cape_mean_annual/DJF/MAM/JJA/SON`: Mean CAPE (J/kg)
  - `potevap_sum_annual/DJF/MAM/JJA/SON`: Total potential evaporation (mm/month)

**Processing rules**:
- Processes 0.125° NLDAS gridded data to county-level aggregates
- Uses spatial weight matrix to map grid cells to counties
- Calculates annual and seasonal (DJF/MAM/JJA/SON) statistics
- Temperature converted from Kelvin to Celsius
- Wind speed calculated from U and V components
- Relative humidity calculated from specific humidity, temperature, and pressure

---

## Conventions

- **Keys**:
  - Use `COUNTY_FIPS` (5-digit string) for county joins.
  - Panel datasets include `Year` (int).
- **Paths**:
  - Original: `Data/Original/...`
  - Processed: `Data/Processed/...`
- **Cleaning assumptions**:
  - Drop all-empty rows at read time.
  - Standardize FIPS to 5-digit string.

---

## Maintenance Notes

After each cleaning script is stabilized, please update this README with:
- Data source description, coverage window, key variables and units.
- Cleaning rules and assumptions.
- Output structure (granularity, row counts if helpful) and join keys to other tables.

## Recent Updates

**Additional Utility Scripts Available**:
- `Code/Clean/CDC_Data_Integrity_Checker.py`: Validates CDC data completeness
- `Code/Clean/PNSP_Unique_Name.py`: Extracts unique pesticide compound names
- `Code/Clean/ENV_GEE_*.py`: Google Earth Engine data collection scripts
- `Code/Test/BSTM_Run.py`: Bayesian spatio-temporal model runner
- `Code/Test/Data_Loading.py`: Model data loading utilities

**File Organization**:
- All cleaning scripts follow consistent path resolution via `config.yaml`
- Environmental data processing includes multiple GEE utility scripts
- Test/modeling scripts are organized in `Code/Test/` directory
- Model results are systematically stored in `Result/Test/` with timestamped summaries

---

## PCA-Derived Covariates (SVI and Climate Factors)

### 1) Master Covariate File

- **Generation script**: `Code/PCA/PCA_Systematic.py`
- **Output file**: `Data/Processed/PCA/Master_Covariates.csv`
- **Granularity**: County × Year (panel; 1999–2020)
- **Description**: This file contains the final, model-ready covariates derived from the systematic PCA workflow. It is the recommended source of covariates for all subsequent analyses, including the Bayesian models.

- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `SVI_PCA`: **Socioeconomic Vulnerability Index**. A composite index where higher values indicate higher socioeconomic vulnerability (e.g., higher poverty, higher unemployment, lower education). Derived from a VIF-screened set of socioeconomic variables.
  - `Climate_Factor_1`, `Climate_Factor_2`, etc.: **Composite Climate Factors**. Orthogonal (uncorrelated) climate indices derived from a VIF-screened set of annual and seasonal meteorological variables. The interpretation of each factor (e.g., "Warm & Dry" vs. "Cool & Wet") can be found in the PCA diagnostic tables in `Result/Tables/PCA_Systematic/`.
  - Other demographic and land use variables included for direct use in modeling.

**Use**: This is the primary source file for covariates in all modeling stages. Join with outcome data via `COUNTY_FIPS` and `Year`.