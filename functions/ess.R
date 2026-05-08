## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/BRAR_gsd
##
## Purpose of script: Functions relating to ESS tables
##
## Author: Corey Voller
##
## Date Created: 17-03-2025
##
## QC'd by:
## QC date:
##
## ─────────────────────────────────────────────────────────────────────────────
##
## Notes:
##   
##
## ─────────────────────────────────────────────────────────────────────────────
##
## 
## 
## Process Expected sample size ------------------------------------------------

# To be used to compare frequentist and bayesian columns
process_EN <- function(EN_group, EN_brar) {
  EN_brar_inf <- EN_brar[, K, , ]  # Extract slice where k = 5
  num_slices <- dim(EN_brar_inf)[3] 
  
  # Compute tidy CI for each slice
  EN_brar_tidy_list <- lapply(1:num_slices, \(z) {
    apply(EN_brar_inf[, , z], 2, \(x) tidy_se(x, digit = dig))
  })
  
  # Compute comparisons for each slice
  comparison_results <- lapply(1:num_slices, \(z) {
    sapply(1:ncol(EN_group), \(j) compare_se(EN_group[, j], EN_brar_inf[, j, z]))
  })
  
  list(tidy = EN_brar_tidy_list, comparison = comparison_results)
}