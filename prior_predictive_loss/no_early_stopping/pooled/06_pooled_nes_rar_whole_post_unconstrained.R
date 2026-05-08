## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/optim-rar-gs
##
## Purpose of script: Losses using whole posterior with no early stopping
##                    with unconstrained allocation
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
# Simulate RAR Pooled ----------------------------------------------------------
simulate_BRAR_pl_thetapr_alloc_negloss <- function(K,
                                                   theta,
                                                   sigma,
                                                   I_theta_fix,
                                                   delta,
                                                   theta_0,
                                                   theta_tau_0) {
  # Draw from prior distribution
  theta_true <- rnorm(1, mean = theta_0, sd = sqrt(1 / theta_tau_0))
  mu_2 <- 0
  # Ratio, estimates of mu, theta
  ratio <- r <- numeric(K)
  # data from trial
  x_1 <- x_2 <- mu1 <- mu2 <- theta_hat <- numeric(K)
  # Information
  I1 <- I2 <- Ip <- numeric(K)
  # Number of patients per arm
  n1 <- n2 <- n1.new <- n2.new <- numeric(K)
  # Theta posterior
  theta_hat_tau<- theta_hat_n <- numeric(K)
  tau_post <- mu_post <- numeric(K)
  # precision
  tau  <- 1 / sigma ^ 2
  # Early stopping vars
  z <- reject <- numeric(K)
  # Inflation factor for GSD
  I_max <- I_theta_fix * inflation_factor
  # Posterior expected loss
  expected_loss <- numeric(K)
  w_exp <- log(a)/delta
  neg_n1 <- 0
  neg_n2 <- 0
  S1 <- S2 <- 0
  # Expected loss at each stage
  loss_stage_post <- loss_stage_real <- numeric(K)
  stop_k <- NA_integer_
  for (k in 1:K) {
    Ip[k] <- (k / K) * I_theta_fix
    if (k == 1) {
      # Initial equal allocation
      ratio[k] <- 1
      r[k] <- 1
      # Enroll 20 patients when N = 100, K=5
      n1.new[k] <- n2.new[k] <- 2 * Ip[k]
      I1[k] <- I2[k] <- 2 * Ip[k]
      dI1 <- I1[k]
      dI2 <- I2[k]
      # simulate initial block means
      x1 <- rnorm(1, theta_true, sqrt(sigma^2 / dI1))
      x2 <- rnorm(1, mu_2,     sqrt(sigma^2 / dI2))
      
      mu1[k] <- x1
      mu2[k] <- x2
      
    } else {
      # unconstrained info allocation (can go up or down)
      I1[k] <- Ip[k] * (1 + r[k])
      I2[k] <- Ip[k] * (1 + (1 / r[k]))
      
      dI1 <- I1[k] - I1[k - 1]
      dI2 <- I2[k] - I2[k - 1]
      
      ## arm 1
      if (dI1 > 0) {
        x1 <- rnorm(1, theta_true, sqrt(sigma^2 / dI1))
        mu1[k] <- (I1[k - 1] * mu1[k - 1] + dI1 * x1) / I1[k]
      } else {
        # negative or zero increment: no new data, mean unchanged
        mu1[k] <- mu1[k - 1]
      }
      
      ## arm 2
      if (dI2 > 0) {
        x2 <- rnorm(1, 0, sqrt(sigma^2 / dI2))
        mu2[k] <- (I2[k - 1] * mu2[k - 1] + dI2 * x2) / I2[k]
      } else {
        mu2[k] <- mu2[k - 1]
      }
    }
    
    theta_hat[k] <- mu1[k] - mu2[k]
    theta_hat_tau[k] <- theta_tau_0 + Ip[k]
    if(theta_hat_tau[k] <= 0) stop("Error: Tau must be positive")
    # Posterior mean
    theta_hat_n[k] <- (theta_tau_0 * theta_0 + theta_hat[k]*Ip[k])/theta_hat_tau[k]
    
    # Posterior expected loss
    mu_post <- theta_hat_n[k]
    posterior_var <- 1/theta_hat_tau[k]
    sigma_post <- sqrt(posterior_var)
    
    p_pos <- pnorm(mu_post / sigma_post)    
    p_neg <- 1 - p_pos                      
    
    # closed-form exponential-weighted integrals
    C2_pos <- exp(w_exp * mu_post + 0.5 * w_exp^2 * posterior_var) *
      pnorm((mu_post + w_exp * posterior_var) / sigma_post)   # B
    
    C1_neg <- exp(-w_exp * mu_post + 0.5 * w_exp^2 * posterior_var) *
      pnorm((w_exp * posterior_var - mu_post) / sigma_post)  # C
    
    C1 <- p_pos + C1_neg
    C2 <- C2_pos + p_neg
    
    loss_stage_post[k] <- I1[k] * p_pos + I2[k] * p_neg +
      I2[k] * C2_pos + I1[k] * C1_neg
    # Realised loss at true theta
    loss_stage_real[k] <- loss_function(
      theta = theta_true,
      I1 = I1[k],
      I2 = I2[k],
      a = a,
      delta = delta
    )
    
    if (k < K) {
      r[k + 1] <- sqrt(C2 / C1)
      # Compute C1 and C2
      ratio[k + 1] = a ** (theta_hat_n[k] / (2 * delta))
    }
    
  }
  if (is.na(stop_k)) stop_k <- K
  c(final_post_loss = loss_stage_post[stop_k],
    final_realised_loss = loss_stage_real[stop_k])
}


run_rar_prior_bayesloss_chunked <- function(
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
          simulate_BRAR_pl_thetapr_alloc_negloss(
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
    "simulate_BRAR_pl_thetapr_alloc_negloss",
    "priors",
    "inflation_factor",
    "I_theta_fix", "delta", "a","a_crit","b_crit"
  ),
  envir = environment()
)

chunk_size <- 100000L
start.time <- Sys.time()

res <- run_rar_prior_bayesloss_chunked(
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

bayestabrarnswpplnegpat <- data.frame(
  prior = 1:nrow(priors),
  theta_0 = priors$theta_0,
  theta_tau0 = priors$theta_tau0,
  est_post = res$est_post,
  mcse_post = res$mcse_post,
  est_real = res$est_real,
  mcse_real = res$mcse_real
)
print(bayestabrarnswpplnegpat)
# saveRDS(bayestabrarnswpplnegpat, file = file.path("output/expected_losses", "bayestabrarnswpplnegpat.rds"))




