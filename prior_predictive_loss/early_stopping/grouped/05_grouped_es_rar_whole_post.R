## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/optim-rar-gs
##
## Purpose of script: Losses for Whole Posterior with Early Stopping
##
## Author: Corey Voller
##
## Date Created: 16-02-2026
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
# Simulate BRAR Grouped early stopping -----------------------------------------
message("Simulate Bayesian RAR Grouped Loss es - Whole Posterior")
simulate_BRAR_grp_thetapr_alloc_final_loss <- function(
    K, theta, sigma, I_theta_fix, delta, theta_0, theta_tau_0
) {
  theta_true <- rnorm(1, mean = theta_0, sd = sqrt(1 / theta_tau_0))
  mu_2 <- 0
  ratio <- r <- numeric(K)
  
  x_1 <- x_2 <- muhat_1 <- muhat_2 <- theta_hat <- theta_hat_group <- numeric(K)
  
  I1 <- I2 <- Ip <- numeric(K)
  
  n1 <- n2 <- numeric(K)
  
  theta_hat_tau <- theta_hat_n <- numeric(K)
  
  # early stopping
  z_stat <- numeric(K)
  stop_k <- NA_integer_
  
  I_max <- I_theta_fix * inflation_factor
  
  w_exp <- log(a) / delta
  # Expected loss at each stage
  loss_stage_post <- loss_stage_real <- numeric(K)
  for (k in 1:K) {
    if (k == 1) {
      ratio[k] <- 1
      r[k] <- 1
      n1[k] <- n2[k] <- (k / K) * 2 * I_max
    } else {
      n1[k] <- (1 / K) * (sigma^2) * I_max * (1 + r[k])
      n2[k] <- (1 / K) * (sigma^2) * I_max * (1 + (1 / r[k]))
      
    }
    # Data draws
    muhat_1[k] <- rnorm(1, mean = theta_true, sd = sqrt(sigma^2 / n1[k])) 
    muhat_2[k] <- rnorm(1, mean = mu_2,     sd = sqrt(sigma^2 / n2[k])) 
    
    I1[k] <- n1[k] / sigma^2
    I2[k] <- n2[k] / sigma^2
    Ip[k] <- 1 / (1 / I1[k] + 1 / I2[k])
    
    theta_hat[k] <- muhat_1[k] - muhat_2[k]
    theta_hat_group[k] = (1 / k) * sum(theta_hat[1:k])
    
    theta_hat_tau[k] <- theta_tau_0 + Ip[k]
    if (theta_hat_tau[k] <= 0) stop("Error: Tau must be positive")
    
    theta_hat_n[k] <- (theta_tau_0 * theta_0 + theta_hat_group[k] * Ip[k]) / theta_hat_tau[k]
    
    # Posterior expected loss 
    mu_post <- theta_hat_n[k]
    posterior_var <- 1 / theta_hat_tau[k]
    sigma_post <- sqrt(posterior_var)
    
    p_pos <- pnorm(mu_post / sigma_post)
    p_neg <- 1 - p_pos
    
    C2_pos <- exp(w_exp * mu_post + 0.5 * w_exp^2 * posterior_var) *
      pnorm((mu_post + w_exp * posterior_var) / sigma_post)
    
    C1_neg <- exp(-w_exp * mu_post + 0.5 * w_exp^2 * posterior_var) *
      pnorm((w_exp * posterior_var - mu_post) / sigma_post)
    
    #  Posterior expected loss
    loss_stage_post[k] <- I1[k] * p_pos + I2[k] * p_neg + I2[k] * C2_pos + I1[k] * C1_neg
    # Realised loss at true theta
    loss_stage_real[k] <- loss_function(
      theta = theta_true,
      I1 = I1[k],
      I2 = I2[k],
      a = a,
      delta = delta
    )
    
    # Early stopping rule
    z_stat[k] <- theta_hat_group[k] * sqrt(sum(Ip[1:k]))
    if (z_stat < a_crit[k] || z_stat > b_crit[k]) { stop_k <- k; break }
    # Allocation update for next stage
    if (k < K) {
      C1 <- p_pos + C1_neg
      C2 <- C2_pos + p_neg
      
      r[k + 1] <- sqrt(C2 / C1)
      ratio[k + 1] <- a ** (theta_hat_n[k] / (2 * delta))
    }
  }
  
  if (is.na(stop_k)) stop_k <- K
  c(final_post_loss = sum(loss_stage_post),
    final_realised_loss = sum(loss_stage_real))
}

run_bayes_prior_loss_chunked <- function(
    cl,
    priors,
    nsims,
    chunk_size = 100000L,
    K = 5L,
    sigma = 1,
    I_theta_fix,
    delta
) {
  Z <- nrow(priors)
  n_chunks <- ceiling(nsims / chunk_size)
  
  est_post <- est_real <- numeric(Z)
  mcse_post <- mcse_real <- numeric(Z)
  
  pb <- txtProgressBar(min = 0, max = Z * n_chunks, style = 3)
  step <- 0L
  
  for (z in 1:Z) {
    pr <- priors[z, ]
    theta_0 <- pr[["theta_0"]]
    theta_tau_0 <- pr[["theta_tau0"]]
    
    S_post <- 0.0; Q_post <- 0.0
    S_real <- 0.0; Q_real <- 0.0
    n_tot <- 0L
    
    for (c in 1:n_chunks) {
      offset <- (c - 1L) * chunk_size
      n_this <- min(chunk_size, nsims - offset)
      
      losses_list <- parLapply(
        cl, 1:n_this,
        function(i, z, offset, K, sigma, I_theta_fix, delta, theta_0, theta_tau_0) {
          set.seed(seed_ij(i, z, offset))
          simulate_BRAR_grp_thetapr_alloc_final_loss(
            K = K,
            sigma = sigma,
            I_theta_fix = I_theta_fix,
            delta = delta,
            theta_0 = theta_0,
            theta_tau_0 = theta_tau_0
          )
        },
        z = z,
        offset = offset,
        K = K, sigma = sigma,
        I_theta_fix = I_theta_fix,
        delta = delta,
        theta_0 = theta_0,
        theta_tau_0 = theta_tau_0
      )
      
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
  list(est_post = est_post,
       mcse_post = mcse_post,
       est_real = est_real,
       mcse_real = mcse_real)
}


num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c(
    "seed_ij", "BASE_SEED", "STRIDE",
    "expected_loss_mu_tau_I",
    "loss_function",
    "simulate_BRAR_grp_thetapr_alloc_final_loss",
    "priors",
    "inflation_factor",
    "I_theta_fix", "delta", "a","a_crit","b_crit"
  ),
  envir = environment()
)

chunk_size <- 100000L
start.time <- Sys.time()

res <- run_bayes_prior_loss_chunked(
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

bayestabrareswpgrp <- data.frame(
  prior = 1:nrow(priors),
  theta_0 = priors$theta_0,
  theta_tau0 = priors$theta_tau0,
  est_post = res$est_post,
  mcse_post = res$mcse_post,
  est_real = res$est_real,
  mcse_real = res$mcse_real
)
print(bayestabrareswpgrp)
# saveRDS(bayestabrareswpgrp, file = file.path("output/expected_losses", "bayestabrareswpgrp.rds"))
