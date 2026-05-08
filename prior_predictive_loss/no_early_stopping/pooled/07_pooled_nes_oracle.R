## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/corey/Documents/PhD/projects/optim-rar-gs
##
## Purpose of script: Oracle No Early stopping
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
# Oracle functions -------------------------------------------------------------
oracle_ratio_abs <- function(theta, a, delta) a^(abs(theta)/(2*delta))
oracle_ratio_dir <- function(theta, a, delta) a^(theta/(2*delta))  # <1 if theta<0

# Simulate RAR Pooled Loss - Oracle NES ----------------------------------------
message("Simulate RAR Pooled Loss - Oracle NES")
simulate_oracle_freqrisk_final_loss <- function(
    K, sigma, I_theta_fix, delta, theta_0, theta_tau_0,
    oracle = c("dir","abs")
) {
  theta_true <- rnorm(1, mean = theta_0, sd = sqrt(1 / theta_tau_0))
  mu_2 <- 0
  ratio <- numeric(K)
  
  oracle <- match.arg(oracle)
  
  r <- if (oracle=="dir") oracle_ratio_dir(theta_true, a, delta) else oracle_ratio_abs(theta_true, a, delta)
  x1 <- x2 <- numeric(K)
  mu1 <- mu2 <- numeric(K)
  theta_hat <- numeric(K)
  tau_post <- mu_post <- numeric(K)
  I1 <- I2 <- Ip <- numeric(K)
  n1 <- n2 <- n1.new <- n2.new <- numeric(K)
  
  I_max <- I_theta_fix * inflation_factor
  w_exp <- log(a) / delta
  # Expected loss at each stage
  loss_stage_post <- loss_stage_real <- numeric(K)
  stop_k <- NA_integer_
  
  for (k in 1:K) {
    if (k == 1) {
      n1.new[k] <- n2.new[k] <- (k / K) * 2 * I_theta_fix
      n1[k] <- n2[k] <- (k / K) * 2 * I_theta_fix
    } else {
      n1[k] <- (k / K) * (sigma^2) * I_theta_fix * (1 + r)
      n2[k] <- (k / K) * (sigma^2) * I_theta_fix * (1 + (1 / r))
      
      n1.new[k] <- n1[k] - n1[k - 1]
      n2.new[k] <- n2[k] - n2[k - 1]
      
      if (n1.new[k] <= 0) {
        n1[k] <- n1[k - 1] + 1e-5
        n1.new[k] <- n1[k] - n1[k - 1]
        n2[k] <- (k * n1[k] * sigma^2 * I_theta_fix) / (5 * n1[k] - k * sigma^2 * I_theta_fix)
        n2.new[k] <- n2[k] - n2[k - 1]
      }
      if (n2.new[k] <= 0) {
        n2[k] <- n2[k - 1] + 1e-5
        n2.new[k] <- n2[k] - n2[k - 1]
        n1[k] <- (k * n2[k] * sigma^2 * I_theta_fix) / (5 * n2[k] - k * sigma^2 * I_theta_fix)
        n1.new[k] <- n1[k] - n1[k - 1]
      }
    }
    
    x1[k] <- rnorm(1, mean=theta_true, sd=sqrt(sigma^2/n1.new[k]))
    x2[k] <- rnorm(1, mean=mu_2,         sd=sqrt(sigma^2/n2.new[k]))
    
    mu1[k] <- sum(n1.new[1:k]*x1[1:k]) / sum(n1.new[1:k])
    mu2[k] <- sum(n2.new[1:k]*x2[1:k]) / sum(n2.new[1:k])
    
    I1[k] <- n1[k] / sigma^2
    I2[k] <- n2[k] / sigma^2
    Ip[k] <- 1 / (1 / I1[k] + 1 / I2[k])
    
    theta_hat[k] <- mu1[k] - mu2[k]
    # Posterior for theta given data up to k
    tau_post[k] <- theta_tau_0 + Ip[k]
    mu_post[k]  <- (theta_tau_0 * theta_0 + Ip[k] * theta_hat[k]) / tau_post[k]
    
    # Posterior expected loss (Bayes loss)
    loss_stage_post[k] <- expected_loss_mu_tau_I(
      mu_post = mu_post[k],
      tau_post = tau_post[k],
      I1 = I1[k],
      I2 = I2[k],
      a = a,
      delta = delta
    )
    loss_stage_real[k] <- loss_function(
      theta = theta_true,
      I1 = I1[k],
      I2 = I2[k],
      a = a,
      delta = delta
    )
  }
  
  if (is.na(stop_k)) stop_k <- K
  c(final_post_loss = loss_stage_post[K],
    final_realised_loss = loss_stage_real[K])
}

# Run 
run_oracle_freqrisk_chunked <- function(
    cl, priors, nsims,
    oracle = c("dir","abs"),
    chunk_size = 100000L,
    K = 5L, sigma = 1,
    I_theta_fix, delta
) {
  oracle <- match.arg(oracle)
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
      offset <- (c-1L)*chunk_size
      n_this <- min(chunk_size, nsims-offset)
      
      losses_list <- parLapply(
        cl, 1:n_this,
        function(i, z, offset, K, sigma, I_theta_fix, delta, oracle, theta_0, theta_tau_0) {
          set.seed(seed_ij(i, z, offset))
          simulate_oracle_freqrisk_final_loss(
            K=K,
            sigma=sigma,
            I_theta_fix=I_theta_fix,
            delta=delta,
            oracle=oracle,
            theta_0 = theta_0,
            theta_tau_0 = theta_tau_0
          )
        },
        z=z, offset=offset,
        K=K, sigma=sigma, I_theta_fix=I_theta_fix, delta=delta, oracle=oracle,
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

num_cores <- max(1, detectCores()-1)
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c(
    "seed_ij","BASE_SEED","STRIDE",
    "oracle_ratio_abs","oracle_ratio_dir",
    "expected_loss_mu_tau_I",
    "loss_function",
    "simulate_oracle_freqrisk_final_loss",
    "theta","nsims","priors",
    "I_theta_fix","delta",
    "a","inflation_factor","a_crit","b_crit","K"
  ),
  envir = environment()
)

chunk_size <- 100000L

res_dir <- run_oracle_freqrisk_chunked(
  cl,
  priors = priors,
  nsims = nsims,
  oracle = "dir",
  chunk_size = chunk_size,
  K = K,
  sigma = 1,
  I_theta_fix = I_theta_fix,
  delta = delta
)

# res_abs <- run_oracle_freqrisk_chunked(
#   cl, theta_vec = theta, nsims = nsims,
#   oracle = "abs",
#   chunk_size = chunk_size,
#   N = 100L, K = K, sigma = 1,
#   I_theta_fix = I_theta_fix, delta = delta
# )

stopCluster(cl)

print(end.time - start.time)

bayestabraroraclenspl <- data.frame(
  prior = 1:nrow(priors),
  theta_0 = priors$theta_0,
  theta_tau0 = priors$theta_tau0,
  est_post = res_dir$est_post,
  mcse_post = res_dir$mcse_post,
  est_real = res_dir$est_real,
  mcse_real = res_dir$mcse_real
)
print(bayestabraroraclenspl)

#print(res_dir$est); print(res_dir$mcse)
#print(res_abs$est); print(res_abs$mcse)
# saveRDS(bayestabraroraclenspl, file = file.path("output/expected_losses", "bayestabraroraclenspl.rds"))