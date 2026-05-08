## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/corey/Documents/PhD/projects/optim-rar-gs
##
## Purpose of script: Frequentist fixed 1:1 With early stopping
##
## Author: Corey Voller
##
## Date Created: 18-12-2025
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
# Simulate fixed 1:1 Grouped ---------------------------------------------------
message("Simulate Fixed 1:1 Group Loss")
simulate_fixed11grp_prior_bayesloss_final <- function(K,
                                                      mu_1,
                                                      mu_2,
                                                      sigma,
                                                      I_theta_fix,
                                                      delta,
                                                      theta_0,
                                                      theta_tau_0) {
  # Draw theta for data generation
  theta_true <- rnorm(1, mean = theta_0, sd = sqrt(1 / theta_tau_0))
  mu_2 <- 0
  # Ratio, estimates of mu, theta
  muhat_1 <- muhat_2 <- theta_hat <- theta_hat_group <- ratio <- numeric(K)
  # Data from the trial (stage means)
  n1 <- n2 <- I1 <- I2 <- numeric(K)
  # posterior objects
  tau_post <- mu_post <- numeric(K)
  # Expected loss at each stage
  loss_stage_post <- loss_stage_real <- numeric(K)
  for (k in 1:K) {
    if (k == 1) {
      ratio[k] <- 1
      # Initial equal allocation
      n1[k] <- n2[k] <- (k / K) * 2 * I_theta_fix
    } else {
      ratio[k] <-  1
      n1[k] <- (1 / K) * (sigma ^ 2) * I_theta_fix * (1 + ratio[k])
      n2[k] <- (1 / K) * (sigma ^ 2) * I_theta_fix * (1 + (1 / ratio[k]))
    }
    # Sample once from normal distribution
    muhat_1[k] <- c(rnorm(1, mean = theta_true, sd = sqrt(sigma ^ 2 / n1[k])))
    muhat_2[k] <- c(rnorm(1, mean = mu_2, sd = sqrt((sigma ^ 2 / n2[k]))))
    # Calculate treatment effect estimate
    theta_hat[k] <- muhat_1[k] - muhat_2[k]
    theta_hat_group[k] = (1 / k) * sum(theta_hat[1:k])
    # Convert to information
    I1[k] <- n1[k] / sigma ^ 2
    I2[k] <- n2[k] / sigma ^ 2
    Ip_k <- 1 / (1 / I1[k] + 1 / I2[k])
    
    # Posterior for theta
    tau_post[k] <- theta_tau_0 + Ip_k
    mu_post[k]  <- (theta_tau_0 * theta_0 + Ip_k * theta_hat[k]) / tau_post[k]
    
    # Posterior expected loss
    loss_stage_post[k] <- expected_loss_mu_tau_I(
      mu_post = mu_post[k],
      tau_post = tau_post[k],
      I1 = I1[k],
      I2 = I2[k],
      a = a,
      delta = delta
    )
    # Realised loss at true theta
    loss_stage_real[k] <- loss_function(
      theta = theta_true,
      I1 = I1[k],
      I2 = I2[k],
      a = a,
      delta = delta
    )
  }
  c(
    final_post_loss = sum(loss_stage_post),
    final_realised_loss = sum(loss_stage_real)
  )
}

# Run function
run_fixed11grp_prior_bayesloss_chunked <- function(cl,
                                                   priors,
                                                   nsims,
                                                   chunk_size = 100000L,
                                                   K = 5L,
                                                   sigma = 1,
                                                   I_theta_fix,
                                                   delta) {
  Z <- nrow(priors)
  n_chunks <- ceiling(nsims / chunk_size)
  
  est_post <- est_real <- numeric(Z)
  mcse_post <- mcse_real <- numeric(Z)
  
  pb <- txtProgressBar(min = 0,
                       max = Z * n_chunks,
                       style = 3)
  step <- 0L
  
  for (z in 1:Z) {
    pr <- priors[z, ]
    theta_0 <- pr[["theta_0"]]
    theta_tau_0 <- pr[["theta_tau0"]]
    
    S_post <- 0.0
    Q_post <- 0.0
    S_real <- 0.0
    Q_real <- 0.0
    n_tot <- 0L
    
    for (c in 1:n_chunks) {
      offset <- (c - 1L) * chunk_size
      n_this <- min(chunk_size, nsims - offset)
      
      losses_list <- parLapply(cl, 1:n_this, function(i,
                                                      z,
                                                      offset,
                                                      K,
                                                      sigma,
                                                      I_theta_fix,
                                                      delta,
                                                      theta_0,
                                                      theta_tau_0) {
        set.seed(seed_ij(i, z, offset))
        simulate_fixed11grp_prior_bayesloss_final(
          K = K,
          sigma = sigma,
          I_theta_fix = I_theta_fix,
          delta = delta,
          theta_0 = theta_0,
          theta_tau_0 = theta_tau_0
        )
      }, z = z, offset = offset, K = K, sigma = sigma, I_theta_fix = I_theta_fix, delta = delta, theta_0 = theta_0, theta_tau_0 = theta_tau_0)
      
      losses_mat <- do.call(rbind, losses_list)
      
      post <- losses_mat[, "final_post_loss"]
      real <- losses_mat[, "final_realised_loss"]
      
      S_post <- S_post + sum(post)
      Q_post <- Q_post + sum(post * post)
      
      S_real <- S_real + sum(real)
      Q_real <- Q_real + sum(real * real)
      
      n_tot <- n_tot + length(post)
      
      step <- step + 1L
      setTxtProgressBar(pb, step)
    }
    
    # Means
    est_post[z] <- S_post / n_tot
    est_real[z] <- S_real / n_tot
    
    # MCSEs
    if (n_tot > 1L) {
      s2_post <- (Q_post - (S_post * S_post) / n_tot) / (n_tot - 1L)
      s2_real <- (Q_real - (S_real * S_real) / n_tot) / (n_tot - 1L)
      mcse_post[z] <- sqrt(s2_post / n_tot)
      mcse_real[z] <- sqrt(s2_real / n_tot)
    } else {
      mcse_post[z] <- NA_real_
      mcse_real[z] <- NA_real_
    }
  }
  
  close(pb)
  list(
    est_post = est_post,
    mcse_post = mcse_post,
    est_real = est_real,
    mcse_real = mcse_real
  )
}

# Parallel coding
num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)
# Export clusters
clusterExport(
  cl,
  varlist = c(
    "seed_ij",
    "BASE_SEED",
    "STRIDE",
    "expected_loss_mu_tau_I",
    "loss_function",
    "simulate_fixed11grp_prior_bayesloss_final",
    "priors",
    "I_theta_fix",
    "delta",
    "a"
  ),
  envir = environment()
)

chunk_size <- 100000L
start.time <- Sys.time()

resgrp <- run_fixed11grp_prior_bayesloss_chunked(
  cl = cl,
  priors = priors,
  nsims = nsims,
  chunk_size = chunk_size,
  K = 5L,
  sigma = 1,
  I_theta_fix = I_theta_fix,
  delta = delta
)

end.time <- Sys.time()
stopCluster(cl)

print(end.time - start.time)

freqtabfixgrp <- data.frame(
  prior = 1:nrow(priors),
  theta_0 = priors$theta_0,
  theta_tau0 = priors$theta_tau0,
  est_post = res$est_post,
  mcse_post = res$mcse_post,
  est_real = res$est_real,
  mcse_real = res$mcse_real
)
print(freqtabfixgrp)
#saveRDS(freqtabfixgrp, file = file.path("output/expected_losses", "freqtabfixgrp.rds"))
