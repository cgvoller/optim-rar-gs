## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/optim-rar-gs
##
## Purpose of script: Main file
##
## Author: Corey Voller
##
## Date Created: 21-02-2025
##
## QC'd by:
## QC date:
##
## ─────────────────────────────────────────────────────────────────────────────
##
## Notes:  powershell -command "tree /f > tree.txt" to output dir tree
##   
##
## ─────────────────────────────────────────────────────────────────────────────
##
## 
## Options----------------------------------------------------------------------
# Remove objects
rm(list=ls())

## Output tree -----------------------------------------------------------------

system('powershell -command "tree /f > tree.txt"')

## Source configurations -------------------------------------------------------

# Load library, source functions and themes
source("config.R")

# Parameters -------------------------------------------------------------------
# Number of digits for outputs
dig <- 1
# Type I error
alpha <- 0.025
# Power = 1 - beta
beta <- 0.1
# Constant used in allocation ratio
a <- 4
# Number of analyses
K <- 5
# Fixing information for fixed sample size to give delta
I_theta_fix <- 50
# Delta (where I_theta_fix is 50)
delta <- (qnorm(1-alpha) + qnorm(1-beta))/sqrt(I_theta_fix)
# treatment difference to be used
theta <- seq(-0.5*delta,2*delta,delta/2)
# Mean of control arm
mu_2 <- 0
# Shared s.d. for experimental and control arm
sigma <- 1
# Number of simulations used in Jennison Paper
JenNsims <- 1000000
# Number of simulations
nsims <- 20000000
# Fractions for output
theta_fractions = as.character(paste(MASS::fractions(theta/delta), "\U1D6FF", sep = ""))

### Stopping bounadries --------------------------------------------------------
message("Create Boundaries")

source("boundaries.R")

# Frequentist Tables -----------------------------------------------------------
message("Frequentist Tables")
### Table one: No early stopping -----------------------------------------------
message("Run no early stopping table")

source("jennisonRes/01_tab_no_early_stopping.R")

### Table two: One group delay--------------------------------------------------
message("Run Group delay table")
source("jennisonRes/02_tab_no_early_stopping_delay.R")

### Table three: Early stopping-------------------------------------------------
message("Run Early stopping table")
source("jennisonRes/03_tab_early_stopping.R")

### Table four: Early stopping and one group delay------------------------------
message("Run Early stopping with group delay table")
source("jennisonRes/04_tab_early_stopping_delay.R")

### Table five: Fixed 3:2 sampling ratio ---------------------------------------
message("Run Fixed 3:2 sampling table")
source("jennisonRes/05_tab_fixed_3_2_sampling.R")

# Bayesian Tables --------------------------------------------------------------
message("Run Bayesian Tables")

### Prior/Historical data ------------------------------------------------------
message("Create prior/historical data scenarios for BRAR")

source("priors.R")


### No Early Stopping ----------------------------------------------------------
#### Pooled --------------------------------------------------------------------
source("prior_predictive_loss/no_early_stopping/pooled/01_pooled_nes_fixed_1_1.R")
source("prior_predictive_loss/no_early_stopping/pooled/02_pooled_nes_fixed_optimal_ratio.R")
source("prior_predictive_loss/no_early_stopping/pooled/03_pooled_nes_rar_mle.R")
source("prior_predictive_loss/no_early_stopping/pooled/04_pooled_nes_rar_post_mean.R")
source("prior_predictive_loss/no_early_stopping/pooled/05_pooled_nes_rar_whole_post.R")
source("prior_predictive_loss/no_early_stopping/pooled/06_pooled_nes_rar_whole_post_unconstrained.R")
source("prior_predictive_loss/no_early_stopping/pooled/07_pooled_nes_oracle.R")

#### Grouped -------------------------------------------------------------------
source("prior_predictive_loss/no_early_stopping/grouped/01_grouped_nes_fixed_1_1.R")
source("prior_predictive_loss/no_early_stopping/grouped/02_grouped_nes_fixed_optimal_ratio.R")
source("prior_predictive_loss/no_early_stopping/grouped/03_grouped_nes_rar_mle.R")
source("prior_predictive_loss/no_early_stopping/grouped/04_grouped_nes_rar_post_mean.R")
source("prior_predictive_loss/no_early_stopping/grouped/05_grouped_nes_rar_whole_post.R")
source("prior_predictive_loss/no_early_stopping/grouped/06_grouped_nes_oracle.R")
### Early Stopping -------------------------------------------------------------
#### Pooled --------------------------------------------------------------------
source("prior_predictive_loss/early_stopping/pooled/01_pooled_es_fixed_1_1.R")
source("prior_predictive_loss/early_stopping/pooled/02_pooled_es_fixed_optimal_ratio.R")
source("prior_predictive_loss/early_stopping/pooled/03_pooled_es_rar_mle.R")
source("prior_predictive_loss/early_stopping/pooled/04_pooled_es_rar_post_mean.R")
source("prior_predictive_loss/early_stopping/pooled/05_pooled_es_rar_whole_post.R")
source("prior_predictive_loss/early_stopping/pooled/06_pooled_es_rar_whole_post_unconstrained.R")
source("prior_predictive_loss/early_stopping/pooled/07_pooled_es_rar_stop_probability.R")
source("prior_predictive_loss/early_stopping/pooled/08_pooled_es_oracle.R")

#### Grouped -------------------------------------------------------------------
source("prior_predictive_loss/early_stopping/grouped/01_grouped_es_fixed_1_1.R")
source("prior_predictive_loss/early_stopping/grouped/02_grouped_es_fixed_optimal_ratio.R")
source("prior_predictive_loss/early_stopping/grouped/03_grouped_es_rar_mle.R")
source("prior_predictive_loss/early_stopping/grouped/04_grouped_es_rar_post_mean.R")
source("prior_predictive_loss/early_stopping/grouped/05_grouped_es_rar_whole_post.R")
source("prior_predictive_loss/early_stopping/grouped/06_grouped_es_rar_stop_probability.R")
source("prior_predictive_loss/early_stopping/grouped/07_grouped_es_oracle.R")

# Print system info ------------------------------------------------------------
Sys.info()
sessionInfo()
Sys.time()