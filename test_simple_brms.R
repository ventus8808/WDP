#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(brms)
})

message("[test] Loading data...")
df <- readr::read_csv("Data/df/EQI_AAMR_Interval.csv", show_col_types = FALSE)

# Filter for C00_C97 (all cancers)
df_sub <- df %>%
  filter(Cancer_Type == "C00_C97", 
         EQI_Period == "0005",
         Time_Period == "2006-2010",
         Lag_Years == 5) %>%
  filter(!is.na(AAMR_lower), !is.na(AAMR_upper), !is.na(Smoking_Rate), !is.na(EQI))

message(sprintf("[test] Filtered to %d rows", nrow(df_sub)))

# Convert variables
df_sub$EQI <- factor(paste0("Q", df_sub$EQI), levels = paste0("Q", 1:5))
df_sub$Smoking_Rate_std <- as.numeric(scale(df_sub$Smoking_Rate))
df_sub$State <- as.factor(df_sub$State)
df_sub$AAMR_mid <- (df_sub$AAMR_lower + df_sub$AAMR_upper) / 2

message("[test] Fitting model...")
fit <- brm(AAMR_mid ~ EQI + Smoking_Rate_std + (1|State), 
           data = df_sub, chains = 1, iter = 500, cores = 1, refresh = 0)
message("[test] Success!")
print(summary(fit))
