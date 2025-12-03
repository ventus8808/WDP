#!/usr/bin/env Rscript
# Quick data check for ridgeline test
# Verifies that C00-C97 data exists for the required scenario

suppressPackageStartupMessages({
  library(data.table)
})

cat("========================================\n")
cat("Data Check for Ridgeline Test\n")
cat("========================================\n\n")

# Find project root
project_root <- normalizePath(".")
data_path <- file.path(project_root,
                       "Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate.csv")

cat("Project root:", project_root, "\n")
cat("Data path:", data_path, "\n\n")

if (!file.exists(data_path)) {
  stop("Data file not found!")
}

cat("Loading data...\n")
dt <- fread(data_path)

cat("Total rows:", nrow(dt), "\n")
cat("Total columns:", ncol(dt), "\n\n")

cat("Column names:\n")
print(names(dt))
cat("\n")

# Check for required columns
req_cols <- c("COUNTY_FIPS", "EQI_Period", "Time_Period", "Lag_Years",
              "Cancer_Type", "AAMR_Lower", "AAMR_Upper", "Smoking_Rate",
              "RUCC", "EQI")
miss_cols <- setdiff(req_cols, names(dt))

if (length(miss_cols) > 0) {
  cat("⚠️  Missing columns:", paste(miss_cols, collapse=", "), "\n\n")
} else {
  cat("✓ All required columns present\n\n")
}

# Check cancer types
cat("Available cancer types:\n")
cancers <- sort(unique(dt$Cancer_Type))
print(cancers)
cat("\n")

if ("C00-C97" %in% cancers) {
  cat("✓ C00-C97 found in data\n\n")
} else {
  cat("⚠️  C00-C97 NOT found in data\n\n")
}

# Check EQI periods
cat("Available EQI periods:\n")
eqi_periods <- sort(unique(dt$EQI_Period))
print(eqi_periods)
cat("\n")

# Check AAMR periods
cat("Available AAMR periods:\n")
aamr_periods <- sort(unique(dt$Time_Period))
print(aamr_periods)
cat("\n")

# Check specific scenario
cat("Checking scenario: C00-C97, 2000-2005 EQI, 2006-2010 AAMR\n")
test_dt <- dt[Cancer_Type == "C00-C97" &
              EQI_Period == "2000-2005" &
              Time_Period == "2006-2010"]

cat("Filtered rows:", nrow(test_dt), "\n")

if (nrow(test_dt) > 0) {
  cat("✓ Data exists for test scenario\n")

  # Add State_FIPS if needed
  if (!"State_FIPS" %in% names(test_dt)) {
    test_dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]
  }

  # Remove NA rows
  test_dt <- test_dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]

  # RUCC restriction
  test_dt <- test_dt[RUCC %in% 1:4 | is.na(RUCC)]

  # Complete cases for model
  test_dt <- test_dt[complete.cases(test_dt[, c("Smoking_Rate", "EQI",
                                                  "AAMR_Lower", "AAMR_Upper",
                                                  "State_FIPS")])]

  cat("After filtering:\n")
  cat("  Rows:", nrow(test_dt), "\n")
  cat("  Counties:", length(unique(test_dt$COUNTY_FIPS)), "\n")
  cat("  States:", length(unique(test_dt$State_FIPS)), "\n\n")

  if (nrow(test_dt) < 50) {
    cat("⚠️  WARNING: Less than 50 observations (n=", nrow(test_dt), ")\n")
  } else {
    cat("✓ Sufficient observations (n=", nrow(test_dt), ")\n")
  }

  cat("\nEQI distribution:\n")
  print(table(test_dt$EQI))

  cat("\nAAMR summary:\n")
  cat("  Lower: ", summary(test_dt$AAMR_Lower), "\n")
  cat("  Upper: ", summary(test_dt$AAMR_Upper), "\n")

  cat("\nSmoking_Rate summary:\n")
  print(summary(test_dt$Smoking_Rate))

} else {
  cat("❌ No data found for test scenario!\n")
}

cat("\n========================================\n")
cat("Data check complete\n")
cat("========================================\n")
