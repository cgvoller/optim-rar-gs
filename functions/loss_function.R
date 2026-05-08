## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/BRAR_gsd
##
## Purpose of script: Loss Function
##
## Author: Corey Voller
##
## Date Created: 27-10-2025
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
## Loss function ---------------------------------------------------------------

loss_function <- function(theta, I1, I2, a, delta) {
  if (theta >= 0) {
    return(I1 + a^(theta / delta) * I2)
  } else {
    return(a^(-theta / delta) * I1 + I2)
  }
}
loss_true_theta <- function(theta, I1, I2, a, delta) {
  w <- log(a) / delta
  if (theta >= 0) I1 + I2 * exp(w * theta) else I2 + I1 * exp(-w * theta)
}
## Posterior expected loss -----------------------------------------------------
expected_loss_mu_tau_I <- function(mu_post, tau_post, I1, I2, a, delta) {
  w_exp <- log(a) / delta
  posterior_var <- 1 / tau_post
  sigma_post <- sqrt(posterior_var)
  p_pos <- pnorm(mu_post / sigma_post)
  p_neg <- 1 - p_pos
  
  term_pos_exp <- I2 * exp(w_exp * mu_post + 0.5 * w_exp^2 * posterior_var) *
    pnorm((mu_post + w_exp * posterior_var) / sigma_post)
  
  term_neg_exp <- I1 * exp(-w_exp * mu_post + 0.5 * w_exp^2 * posterior_var) *
    pnorm((w_exp * posterior_var - mu_post) / sigma_post)
  
  I1 * p_pos + I2 * p_neg + term_pos_exp + term_neg_exp
}
