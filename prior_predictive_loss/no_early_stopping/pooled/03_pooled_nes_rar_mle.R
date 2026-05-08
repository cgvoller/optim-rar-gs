## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/optim-rar-gs
##
## Purpose of script: Losses for No early stopping table one (freq pooled+group)
##
## Author: Corey Voller
##
## Date Created: 25-02-2025
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
message("Simulate RAR No Early Stopping Pooled Loss")
simulate_RAR <- function(K,
                         mu_1,
                         mu_2,
                         sigma,
                         I_theta_fix,
                         delta,
                         theta_0,
                         theta_tau_0) {
  # Ratio, estimates of mu, theta
  # Draw theta for data generation
  theta_true <- rnorm(1, mean = theta_0, sd = sqrt(1 / theta_tau_0))
  mu_2 <- 0
  
  muhat_1 <- muhat_2 <- ratio <- numeric(K)
  theta_hat <- rep(NA, K)
  # data from trial
  x_1 <- x_2 <- numeric(0)
  # Information
  I1 <- I2 <- numeric(K)
  # Number of patients per arm
  n1 <- n2 <- n1.new <- n2.new <- numeric(K)
  
  tau_post <- mu_post <- numeric(K)
  # Expected loss at each stage
  loss_stage_post <- loss_stage_real <- numeric(K)
  w_exp <- log(a) / delta
  neg_n1 <- 0
  neg_n2 <- 0
  for (k in 1:K) {
    if (k == 1) {
      # Initial equal allocation
      ratio[k] <- 1
      # Enroll 20 patients when N = 100, K=5
      n1.new[k] <- n2.new[k] <- (k / K) * 2 * I_theta_fix
      n1[k] <- n2[k] <- (k / K) * 2 * I_theta_fix
    } else {
      # Define ratio based on formula with previous estimate of theta
      ratio[k] <-  a ** (theta_hat[k - 1] / (2 * delta))
      # Patients in N1
      n1[k] <- (k / K) * (sigma ^ 2) * I_theta_fix * (1 + ratio[k])
      # Patients in N2
      n2[k] <- (k / K) * (sigma ^ 2) * I_theta_fix * (1 + (1 / ratio[k]))
      # Difference in patients from previous
      n1.new[k] <- n1[k] - n1[k - 1]
      n2.new[k] <- n2[k] - n2[k - 1]
      if (n1.new[k] <= 0) {
        neg_n1 <- neg_n1 + 1
        n1[k] = n1[k - 1] + 0.00001
        n1.new[k] = n1[k] - n1[k - 1]
        n2[k] = (k * n1[k] * sigma ^ 2 * I_theta_fix) / (5 * n1[k] - k *
                                                           sigma ^ 2 * I_theta_fix)
        n2.new[k] = n2[k] - n2[k - 1]
      }
      if (n2.new[k] <= 0) {
        neg_n2 <- neg_n2 + 1
        n2[k] = n2[k - 1] + 0.00001
        n2.new[k] = n2[k] - n2[k - 1]
        n1[k] = (k * n2[k] * sigma ^ 2 * I_theta_fix) / (5 * n2[k] - k *
                                                           sigma ^ 2 * I_theta_fix)
        n1.new[k] = n1[k] - n1[k - 1]
      }
    }
    # Sample once from normal distribution
    x_1[k] <- rnorm(1, theta_true, sqrt(sigma ^ 2 / n1.new[k]))
    x_2[k] <- rnorm(1, mu_2, sqrt(sigma ^ 2 / n2.new[k]))
    # Calculate treatment effect estimate
    muhat_1[k] <- sum(n1.new[1:k] * x_1[1:k]) / sum(n1.new[1:k])
    muhat_2[k] <- sum(n2.new[1:k] * x_2[1:k]) / sum(n2.new[1:k])
    theta_hat[k] <- muhat_1[k] - muhat_2[k]
    # Convert to information
    I1[k] <- n1[k] / sigma ^ 2
    I2[k] <- n2[k] / sigma ^ 2
    Ip_k <- 1 / (1 / I1[k] + 1 / I2[k])
    
    # Posterior for theta given data up to k
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
  # Return both at final analysis
  c(
    setNames(ratio, paste0("ratio_", 1:K)),
    final_post_loss = loss_stage_post[K],
    final_realised_loss = loss_stage_real[K]
  )
}
run_RAR_prior_bayesloss_chunked <- function(cl,
                                            priors,
                                            nsims,
                                            chunk_size = 100000L,
                                            K = 5L,
                                            sigma = 1,
                                            I_theta_fix,
                                            delta) {
  Z <- nrow(priors)
  n_chunks <- ceiling(nsims / chunk_size)
  ratio_cols <- paste0("ratio_", 1:K)
  
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
        simulate_RAR(
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
      
      keep_n <- max(1L, as.integer(ceiling(ratio_sample_frac * n_this)))
      keep_idx <- if (ratio_sample_frac >= 1)
        seq_len(n_this)
      else
        sample.int(n_this, keep_n)
      
      ratio_dt <- data.table::as.data.table(losses_mat[keep_idx, ratio_cols, drop = FALSE])
      ratio_dt[, prior := z]
      ratio_dt[, sim_id := offset + keep_idx]
      
      # saveRDS(
      #   ratio_dt,
      #   file = file.path(ratio_out_dir, sprintf("prior_%03d_chunk_%05d.rds", z, c))
      # )
      
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


num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c(
    "seed_ij",
    "BASE_SEED",
    "STRIDE",
    "expected_loss_mu_tau_I",
    "loss_function",
    "simulate_RAR",
    "priors",
    "I_theta_fix",
    "delta",
    "a"
  ),
  envir = environment()
)

chunk_size <- 100000L
start.time <- Sys.time()

res <- run_RAR_prior_bayesloss_chunked(
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

freqtabrarnspl <- data.frame(
  prior = 1:nrow(priors),
  theta_0 = priors$theta_0,
  theta_tau0 = priors$theta_tau0,
  est_post = res$est_post,
  mcse_post = res$mcse_post,
  est_real = res$est_real,
  mcse_real = res$mcse_real
)
print(freqtabrarnspl)
# saveRDS(freqtabrarnspl, file = file.path("output/expected_losses", "freqtabrarnspl.rds"))
