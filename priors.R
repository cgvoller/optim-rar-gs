## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/optim-rar-gs
##
## Purpose of script: Generate Priors
##
## Author: Corey Voller 
##
## Date Created: 11-03-2025
##
## QC'd by:
## QC date:
##
## ─────────────────────────────────────────────────────────────────────────────
##
## Notes: A data frame containing the prior values
##        Can extend to a grid of values
##  
## ─────────────────────────────────────────────────────────────────────────────
##
## 
## Generate Priors -------------------------------------------------------------

make_prior_95_0_delta <- function(delta, mass = 0.95) {
  z <- qnorm((1 + mass)/2)         
  sd_theta  <- (delta/2) / z
  tau_theta <- 1 / sd_theta^2
  

  tau01 <- tau02 <- 2 * tau_theta
  
  data.frame(
    mu_1   = delta/2,
    mu_2   = 0,
    tau_01 = tau01,
    tau_02 = tau02
  )
}

new_prior <- make_prior_95_0_delta(delta)

priors <- data.frame(
  mu_1 = c(0, delta, delta, 2 * delta * ((log(3/2)) / log(4)),delta/2,delta/2,delta/2,delta/2),
  mu_2 = rep(0, 8),
  tau_01 = c(0.0000001, 20, 100, 100000,5, 10,20,5),
  tau_02 = c(0.0000001, 20, 100, 100000,5, 10,20,15))

priors <- rbind(priors, new_prior)
priors$theta_tau0 <- 1/(1/priors$tau_01 + 1/priors$tau_02)
priors$theta_0 <- priors$mu_1 - priors$mu_2

rownames(priors) <- paste0("prior ", seq_len(nrow(priors)))

# Filter priors wanted
#priors <- priors[c(5,6,7,9),]
priors <- priors[c(7,9),]