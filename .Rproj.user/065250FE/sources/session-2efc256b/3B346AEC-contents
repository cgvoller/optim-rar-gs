## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/optim-rar-gs
##
## Purpose of script: Create Pampallona and Tsiatis boundaries
##
## Author: Corey Voller
##
## Date Created: 04-04-2025
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
## Boundary --------------------------------------------------------------------
# Design params
design <- rpact::getDesignGroupSequential(
  sided = 1, 
  alpha = 0.025, 
  beta = 0.1,
  typeOfDesign = "PT",       
  deltaPT1 = 0,              
  deltaPT0 = 0,              
  kMax = 5,                 
  #informationRates = c(0.2, 0.4, 0.6, 0.8, 1),
  bindingFutility = TRUE
)
# Futility boundary

a_crit <- c(design$futilityBounds, tail(design$criticalValues, n = 1))

# Efficacy boundary

b_crit <- design$criticalValues

# Inflation factor - to be applied to information levels when using early 
# stopping

inflation_factor <- getDesignCharacteristics(design)$inflationFactor