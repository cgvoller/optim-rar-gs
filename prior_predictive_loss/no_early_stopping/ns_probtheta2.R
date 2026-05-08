## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/corey/Documents/PhD/projects/BRAR_gsd
##
## Purpose of script: Account for probability of stopping
##
## Author: Corey Voller
##
## Date Created: 06-02-2026
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
############################################################
## STOP-WEIGHTED GOLDEN-SECTION BRAR
##
## Key change vs simulate_BRAR_pl_thetapr_final_loss:
## - ratio[k] (for k>=2) is chosen by golden-section search to minimize
##     int pi(theta | data_{k-1}) P(S_k | data_{k-1}; theta) L(theta; I1^k, I2^k) dtheta
## - BUT loss_stage[k] is still posterior expected loss (closed form),
##   so the final numbers are comparable to previous methods
############################################################
if(RUN_POOLED){
  message("Simulate Whole posterior, accounting for probability of theta + No Early Stopping Pooled Loss")
loss_theta_I <- function(theta, I1, I2, w_exp) {
  ifelse(theta > 0,
         I1 + I2 * exp(w_exp * theta),
         I2 + I1 * exp(-w_exp * theta))
}

############################
## 3) P(S_k | data_{k-1}; theta) via Normal tails on Z_k
## Z_k = theta_hat_k * sqrt(Ip_k)
## Conditional on Z_{k-1}=z_prev and theta:
##   Z_k ~ N( m(theta), s^2 )
## where:
##   dI = Ip_k - Ip_{k-1}
##   m = ( z_prev*sqrt(I_prev) + theta*dI ) / sqrt(I_k)
##   s^2 = dI / I_k
############################
p_stop_at_k_given_theta <- function(theta, z_prev, I_prev, I_k, a_k, b_k) {
  dI <- I_k - I_prev
  if (dI <= 0) return(0)
  
  m <- (z_prev * sqrt(I_prev) + theta * dI) / sqrt(I_k)
  s <- sqrt(dI / I_k)
  
  pnorm((a_k - m) / s) + (1 - pnorm((b_k - m) / s))
}

############################
## 4) Monte Carlo approximation integral at stage k:
##   E_{theta|data_{k-1}}[ P(S_k|theta) * L(theta; I1^k,I2^k) ]
##
## Uses common random numbers (theta_draws fixed) for optimize().
############################
stop_weighted_integral_mc <- function(theta_draws,
                                      z_prev, I_prev,
                                      I_k, a_k, b_k,
                                      I1_k, I2_k,
                                      w_exp) {
  pS <- p_stop_at_k_given_theta(theta_draws, z_prev, I_prev, I_k, a_k, b_k)
  L  <- loss_theta_I(theta_draws, I1_k, I2_k, w_exp)
  mean(pS * L)
}

############################
## 5) Golden-section choice of ratio[k] at stage k (k>=2)
##
## - posterior at k-1: theta | data_{k-1} ~ N(mu_prev, 1/tau_prev)
## - z_prev, I_prev from stage k-1
## - for each candidate ratio r, you set stage-k I1/I2 using YOUR schedule
## - objective is the stop-weighted integral above
choose_ratio_golden <- function(mu_prev, tau_prev,
                                z_prev, I_prev,
                                k, K, I_theta_fix, sigma,
                                a_crit, b_crit,
                                a, delta,
                                M_theta = 4000L,
                                r_bounds = c(1e-3, 1e3),
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  theta_draws <- rnorm(M_theta, mean = mu_prev, sd = 1 / sqrt(tau_prev))
  
  w_exp <- log(a) / delta
  a_k <- a_crit[k]
  b_k <- b_crit[k]
  
  f_logr <- function(logr) {
    r <- exp(logr)
    
    n1_k <- (k / K) * (sigma^2) * I_theta_fix * (1 + r)
    n2_k <- (k / K) * (sigma^2) * I_theta_fix * (1 + (1 / r))
    
    I1_k <- n1_k / sigma^2
    I2_k <- n2_k / sigma^2
    I_k  <- 1 / (1 / I1_k + 1 / I2_k)  # Ip_k
    
    stop_weighted_integral_mc(theta_draws,
                              z_prev, I_prev,
                              I_k, a_k, b_k,
                              I1_k, I2_k,
                              w_exp)
  }
  
  opt <- optimize(f_logr, interval = log(r_bounds))
  exp(opt$minimum)
}

simulate_BRAR_pl_thetapr_final_loss_golden <- function(
    N, K, sigma, I_theta_fix, delta, theta_0, theta_tau_0,
    M_theta = 4000L, r_bounds = c(1e-3, 1e3)
) {
  w_exp <- log(a) / delta
  theta_true <- rnorm(1, mean = theta_0, sd = sqrt(1 / theta_tau_0))
  mu_1 <- theta_true
  mu_2 <- 0
  
  ratio <- numeric(K)
  
  x_1 <- x_2 <- muhat_1 <- muhat_2 <- theta_hat <- numeric(K)
  I1 <- I2 <- Ip <- numeric(K)
  
  n1 <- n2 <- n1.new <- n2.new <- numeric(K)
  
  tau_post <- mu_post <- numeric(K)
  z_stat <- numeric(K)
  loss_stage_post <- loss_stage_real <- numeric(K)
  stop_k <- NA_integer_
  
  for (k in 1:K) {
    if (k == 1) {
      ratio[k] <- 1
      n1.new[k] <- n2.new[k] <- (k / K) * 2 * I_theta_fix
      n1[k] <- n2[k] <- (k / K) * 2 * I_theta_fix
    } else {
      # NEW: choose ratio[k] by golden-section to minimize integral at stage k
      ratio[k] <- choose_ratio_golden(
        mu_prev = mu_post[k - 1],
        tau_prev = tau_post[k - 1],
        z_prev = z_stat[k - 1],
        I_prev = Ip[k - 1],
        k = k, K = K,
        I_theta_fix = I_theta_fix,
        sigma = sigma,
        a_crit = a_crit, b_crit = b_crit,
        a = a, delta = delta,
        M_theta = M_theta,
        r_bounds = r_bounds,
        seed = NULL
      )
      
      n1[k] <- (k / 5) * (sigma^2) * I_theta_fix * (1 + ratio[k])
      n2[k] <- (k / 5) * (sigma^2) * I_theta_fix * (1 + (1 / ratio[k]))
      
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
    
    # Data draws
    x_1[k] <- if (n1.new[k] > 0) rnorm(1, mean = mu_1, sd = sqrt(sigma^2 / n1.new[k])) else NA_real_
    x_2[k] <- if (n2.new[k] > 0) rnorm(1, mean = mu_2, sd = sqrt(sigma^2 / n2.new[k])) else NA_real_
    
    # Cumulative means
    muhat_1[k] <- sum(n1.new[1:k] * x_1[1:k], na.rm = TRUE) / sum(n1.new[1:k])
    muhat_2[k] <- sum(n2.new[1:k] * x_2[1:k], na.rm = TRUE) / sum(n2.new[1:k])
    
    I1[k] <- n1[k] / sigma^2
    I2[k] <- n2[k] / sigma^2
    Ip[k] <- 1 / (1 / I1[k] + 1 / I2[k])
    
    theta_hat[k] <- muhat_1[k] - muhat_2[k]
    
    # Posterior for theta given data up to k
    tau_post[k] <- theta_tau_0 + Ip[k]
    mu_post[k]  <- (theta_tau_0 * theta_0 + Ip[k] * theta_hat[k]) / tau_post[k]
    
    # Posterior expected loss 
    loss_stage_post[k] <- expected_loss_mu_tau_I(
      mu_post = mu_post[k],
      tau_post = tau_post[k],
      I1 = I1[k],
      I2 = I2[k],
      a = a,
      delta = delta
    )
    loss_stage_real[k] <- loss_theta_I(
      theta = theta_true,
      I1 = I1[k],
      I2 = I2[k],
      w_exp
    )
    # Early stopping rule
    z_stat[k] <- theta_hat[k] * sqrt(Ip[k])
    #if (z_stat[k] < a_crit[k]) { stop_k <- k; break }
    #if (z_stat[k] > b_crit[k]) { stop_k <- k; break }
  }
  
  if (is.na(stop_k)) stop_k <- K
  c(final_post_loss = loss_stage_post[stop_k],
    final_realised_loss = loss_stage_real[stop_k])
}

run_bayes_prior_loss_chunked_golden <- function(
    cl,
    priors,
    nsims,
    chunk_size = 10000L,
    N = 100L,
    K = 5L,
    sigma = 1,
    I_theta_fix,
    delta,
    M_theta = 4000L,
    r_bounds = c(1e-3, 1e3)
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
        function(i, z, offset, N, K, sigma, I_theta_fix, delta, theta_0, theta_tau_0, M_theta, r_bounds) {
          set.seed(seed_ij(i, z, offset))
          simulate_BRAR_pl_thetapr_final_loss_golden(
            N = N, K = K,
            sigma = sigma,
            I_theta_fix = I_theta_fix,
            delta = delta,
            theta_0 = theta_0,
            theta_tau_0 = theta_tau_0,
            M_theta = M_theta,
            r_bounds = r_bounds
          )
        },
        z = z,
        offset = offset,
        N = N, K = K, sigma = sigma,
        I_theta_fix = I_theta_fix,
        delta = delta,
        theta_0 = theta_0,
        theta_tau_0 = theta_tau_0,
        M_theta = M_theta,
        r_bounds = r_bounds
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

num_cores <- max(1L, detectCores() - 1L)
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c(
    "seed_ij","BASE_SEED","STRIDE",
    "expected_loss_mu_tau_I",
    "loss_theta_I",
    "p_stop_at_k_given_theta",
    "stop_weighted_integral_mc",
    "choose_ratio_golden",
    "simulate_BRAR_pl_thetapr_final_loss_golden",
    "run_bayes_prior_loss_chunked_golden",
    "priors",
    "inflation_factor",
    "I_theta_fix","delta","a","a_crit","b_crit"
  ),
  envir = environment()
)

chunk_size <- if (exists("chunk_size")) chunk_size else 100000L
start.time <- Sys.time()

res <- run_bayes_prior_loss_chunked_golden(
  cl = cl,
  priors = priors,
  nsims = nsims,
  chunk_size = chunk_size,
  N = 100L,
  K = 5L,
  sigma = 1,
  I_theta_fix = I_theta_fix,
  delta = delta,
  M_theta = 2000L,          # start small, increase for stability
  r_bounds = c(1e-3, 1e3)
)

end.time <- Sys.time()
stopCluster(cl)

print(end.time - start.time)

bayestabrarprobthetanspl <- data.frame(
  prior = 1:nrow(priors),
  theta_0 = priors$theta_0,
  theta_tau0 = priors$theta_tau0,
  est_post = res$est_post,
  mcse_post = res$mcse_post,
  est_real = res$est_real,
  mcse_real = res$mcse_real
)
print(bayestabrarprobthetanspl)

saveRDS(bayestabrarprobthetanspl, file = file.path("output/expected_losses", "bayestabrarprobthetanspl.rds"))
}
