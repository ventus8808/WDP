#!/usr/bin/env Rscript
#'
#' C00-C97 Complete Analysis Runner - Simplified
#'

message("[C00-C97] Starting complete analysis for All Cancers (C00-C97)")

# Key scenarios to run for C00-C97
key_scenarios <- c(
  "AllCancer_TotalEQI_Lag5",              # 1. 总EQI
  "AllCancer_AirEQI_Lag5",                # 2. 空气EQI  
  "AllCancer_WaterEQI_Lag5",              # 3. 水EQI
  "AllCancer_LandEQI_Lag5",               # 4. 土地EQI
  "AllCancer_BuiltEQI_Lag5",              # 5. 建成环境EQI
  "AllCancer_SocioEQI_Lag5",              # 6. 社会人口学EQI
  "AllCancer_RUCC_Urban_TotalEQI_Lag5",   # 7. 城市地区 - 总EQI
  "AllCancer_RUCC_Urban_AirEQI_Lag5",     # 8. 城市地区 - 空气EQI
  "AllCancer_RUCC_Urban_WaterEQI_Lag5",   # 9. 城市地区 - 水EQI
  "AllCancer_RUCC_Rural_TotalEQI_Lag5",   # 10. 农村地区 - 总EQI
  "AllCancer_RUCC_Rural_AirEQI_Lag5"      # 11. 农村地区 - 空气EQI
)

message(sprintf("[C00-C97] Running %d key scenarios", length(key_scenarios)))

success_count <- 0
failure_count <- 0

for (i in seq_along(key_scenarios)) {
  scenario_name <- key_scenarios[i]
  
  message(sprintf("[C00-C97] [%d/%d] Running: %s", i, length(key_scenarios), scenario_name))
  
  # Build and run command
  result <- system2("Rscript", 
                   args = c("Code/brms/02_run_brms_model.R", "--scenario", scenario_name),
                   stdout = FALSE, stderr = FALSE)
  
  if (result == 0) {
    success_count <- success_count + 1
    message(sprintf("[C00-C97] ✓ Completed: %s", scenario_name))
  } else {
    failure_count <- failure_count + 1
    message(sprintf("[C00-C97] ✗ Failed: %s", scenario_name))
  }
}

message(sprintf("[C00-C97] Analysis complete! Success: %d, Failures: %d", 
                success_count, failure_count))

# Generate final report
message("[C00-C97] Generating final report...")
post_result <- system2("Rscript", args = "Code/brms/04_process_results.R", 
                      stdout = FALSE, stderr = FALSE)

if (post_result == 0) {
  message("[C00-C97] ✓ Final report generated successfully")
  message("[C00-C97] Check Results:")
  message("[C00-C97]   - Combined results: Result/brms/reports/brms_fixed_effects_combined.csv")
  message("[C00-C97]   - LMM compatible: Result/brms/brms_lmm_compatible.csv") 
  message("[C00-C97]   - Forest plot: Result/brms/figures/brms_fixed_effects_forest.png")
} else {
  message("[C00-C97] ✗ Report generation failed")
}