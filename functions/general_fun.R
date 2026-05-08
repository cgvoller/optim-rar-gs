## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/BRAR_gsd
##
## Purpose of script: General functions
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
## Notes:
##   
##
## ─────────────────────────────────────────────────────────────────────────────
##
## Seed ------------------------------------------------------------------------
BASE_SEED <- 42L
STRIDE <- 100000000L
seed_ij <- function(i, j, offset = 0L) {
  BASE_SEED + as.integer(offset + i) + as.integer(STRIDE) * as.integer(j - 1L)
}
## Round -----------------------------------------------------------------------
round2 = function(x, digits) {
  posneg = sign(x)
  z = abs(x)*10^digits
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^digits
  z*posneg
}
# Confidence Interval ----------------------------------------------------------
# Create function for CI's of expected sample sizes from simulations
ci_funct <- function(x,digit=3) {
  upper = mean(x) + (qnorm(1 - alpha) * (sd(x) / sqrt(length(x))))
  lower = mean(x) - (qnorm(1 - alpha) * (sd(x) / sqrt(length(x))))
  return(paste0(
    round2(mean(x), digit),
    " [",
    round2(lower, digit),
    ",",
    round2(upper, digit),
    "]"
  ))
}

# Confidence Interval 2 --------------------------------------------------------

# To be used when requiring the upper and lower for comparisons
# e.g., comparing bayesian and frequentist intervals
ci_funct2 <- function(x, alpha = 0.05, digit = 3) {
  upper = mean(x) + (qnorm(1 - alpha) * (sd(x) / sqrt(length(x))))
  lower = mean(x) - (qnorm(1 - alpha) * (sd(x) / sqrt(length(x))))
  return(c(lower = lower, upper = upper, mean = mean(x)))
}

# Tidy CI ----------------------------------------------------------------------

# Formats the CI's for tidy output in tables
tidy_ci <- function(x,digit){
  paste0(
    round2(x["mean"], digit),
    " [",
    round2(x["lower"], digit),
    ",",
    round2(x["upper"], digit),
    "]"
  )
}

# Compare CI's -----------------------------------------------------------------

# Compare two CI's, if their intervals overlap assign equal
compare_cis <- function(ci1, ci2) {
  overlap <- !(ci1[2] < ci2[1] | ci1[1] > ci2[2]) 
  if (overlap) {
    return("equals")
  } else {
    if (ci1[1] > ci2[2]) {
      return("arrow-down")
    } else {
      return("arrow-up")
    }
  }
}

# Standard error ---------------------------------------------------------------

se_funct <- function(x,digit) {
  error <- sd(x) / sqrt(length(x))
  upper = mean(x) + error
  lower = mean(x) - error
  return(paste0(round2(mean(x), digit), " (", round2(error,digit), ")"))
}


se_funct2 <- function(x) {
  error <- sd(x) / sqrt(length(x))
  upper = mean(x) + error
  lower = mean(x) - error
  return(c(lower = lower, upper = upper, mean = mean(x),error=error))
}

tidy_se <- function(x, digit) {
  paste0(round2(x["mean"], digit), " (", round2(x["error"], digit), ")")
}


compare_se <- function(se1, se2) {
  se1 <- round2(se1,dig)
  se2 <- round2(se2,dig)
  overlap <- !(se1[2] < se2[1] | se1[1] > se2[2]) 
  if (overlap) {
    return("equals")
  } else {
    if (se1[1] > se2[2]) {
      return("arrow-down")
    } else {
      return("arrow-up")
    }
  }
}

# Progress bar -----------------------------------------------------------------

progress_fun <- function() {
  progress_env$counter <- progress_env$counter + 1
  setTxtProgressBar(pb, progress_env$counter)
  NULL
}


# Save function ----------------------------------------------------------------
save_outputs <- function(table, name,dir) {
  saveRDS(table, file.path(dir, paste0(name, ".rds")))
  gtsave(table, file.path(dir, paste0(name, ".png")))
  gtsave(table, file.path(dir, paste0(name, ".tex")))
}

# Allocation ratios ------------------------------------------------------------
# Allocation ratio 1 -----------------------------------------------------------

# alloc_ratio <- function(a, theta_hat, delta) {
#   return( a ** (theta_hat/ (2 * delta)))
# }



