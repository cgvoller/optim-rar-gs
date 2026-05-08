## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/corey/Documents/PhD/projects/optim-rar-gs
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

LOSS_MODE <- "terminal"

## Stop-Weighted Golden-Section BRAR -------------------------------------------
required_objs <- c(
  "a",
  "delta",
  "inflation_factor",
  "a_crit",
  "b_crit",
  "priors",
  "I_theta_fix",
  "K",
  "sigma",
  "nsims",
  "seed_ij",
  "BASE_SEED",
  "STRIDE"
)

loss_function <- function(theta, I1, I2, a, delta) {
  if (theta >= 0) {
    I1 + a^(theta / delta) * I2
  } else {
    a^(-theta / delta) * I1 + I2
  }
}

expected_loss_mu_tau_I <- function(mu_post, tau_post, I1, I2, a, delta) {
  w_exp <- log(a) / delta
  v <- 1 / tau_post
  sd_post <- sqrt(v)
  
  p_pos <- pnorm(mu_post / sd_post)
  p_neg <- 1 - p_pos
  
  term_pos_exp <- I2 * exp(w_exp * mu_post + 0.5 * w_exp^2 * v) *
    pnorm((mu_post + w_exp * v) / sd_post)
  
  term_neg_exp <- I1 * exp(-w_exp * mu_post + 0.5 * w_exp^2 * v) *
    pnorm((w_exp * v - mu_post) / sd_post)
  
  I1 * p_pos + I2 * p_neg + term_pos_exp + term_neg_exp
}

r_next_wp <- function(mu_post, tau_post, a, delta, r_bounds = c(0.05, 20)) {
  w <- log(a) / delta
  v <- 1 / tau_post
  sd <- sqrt(v)
  
  p_pos <- pnorm(mu_post / sd)
  p_neg <- 1 - p_pos
  
  C2_pos <- exp(w * mu_post + 0.5 * w^2 * v) *
    pnorm((mu_post + w * v) / sd)
  
  C1_neg <- exp(-w * mu_post + 0.5 * w^2 * v) *
    pnorm((w * v - mu_post) / sd)
  
  C1 <- p_pos + C1_neg
  C2 <- C2_pos + p_neg
  
  r <- if (!is.finite(C1) || !is.finite(C2) || C1 <= 0 || C2 <= 0) 1 else sqrt(C2 / C1)
  max(r_bounds[1], min(r_bounds[2], r))
}

# Probability of stopping ------------------------------------------------------
p_stop_next_given_theta <- function(theta, z_k, I_k, I_next, a_next, b_next) {
  dI <- I_next - I_k
  if (dI <= 0) return(0)
  
  m <- (z_k * sqrt(I_k) + theta * dI) / sqrt(I_next)
  s <- sqrt(dI / I_next)
  
  pnorm((a_next - m) / s) + (1 - pnorm((b_next - m) / s))
}

simulate_BRAR_grp_thetapr_final_loss_golden <- function(
    K, sigma, I_theta_fix, delta, theta_0, theta_tau_0,
    loss_mode = c("accum","terminal"),
    r_bounds = c(0.05, 20),
    return_trace = FALSE
) {
  loss_mode <- match.arg(loss_mode)
  
  # prior predictive true theta
  theta_true <- rnorm(1, mean = theta_0, sd = sqrt(1 / theta_tau_0))
  mu_2 <- 0
  
  I_max <- I_theta_fix * inflation_factor
  Ip_stage <- I_max / K  # fixed pooled increment per stage
  
  # storage
  r <- numeric(K)
  r[1] <- 1
  
  theta_hat_stage <- numeric(K)
  theta_hat_group <- numeric(K)
  
  I1_stage <- I2_stage <- Ip <- numeric(K)
  
  mu_post <- tau_post <- numeric(K)
  z_stat <- numeric(K)
  
  loss_post_stage <- loss_real_stage <- numeric(K)
  stop_k <- NA_integer_
  
  trace <- NULL
  if (return_trace) {
    trace <- tibble()
  }
  
  for (k in 1:K) {
    
    I1_stage[k] <- Ip_stage * (1 + r[k])
    I2_stage[k] <- Ip_stage * (1 + 1 / r[k])
    Ip[k] <- 1 / (1 / I1_stage[k] + 1 / I2_stage[k]) 
    
    n1 <- sigma^2 * I1_stage[k]
    n2 <- sigma^2 * I2_stage[k]
    
    x1 <- rnorm(1, mean = theta_true, sd = sqrt(sigma^2 / n1))
    x2 <- rnorm(1, mean = mu_2,      sd = sqrt(sigma^2 / n2))
    
    theta_hat_stage[k] <- x1 - x2
    theta_hat_group[k] <- mean(theta_hat_stage[1:k])
    
    I_cum <- k * Ip_stage
    
    tau_post[k] <- theta_tau_0 + Ip_stage
    mu_post[k]  <- (theta_tau_0 * theta_0 + Ip_stage * theta_hat_stage[k]) / tau_post[k]
    
    loss_post_stage[k] <- expected_loss_mu_tau_I(mu_post[k], tau_post[k], I1_stage[k], I2_stage[k], a, delta)
    loss_real_stage[k] <- loss_function(theta_true, I1_stage[k], I2_stage[k], a, delta)
    
    z_stat[k] <- theta_hat_group[k] * sqrt(I_cum)
    if (z_stat[k] < a_crit[k] || z_stat[k] > b_crit[k]) { stop_k <- k }
    
    # optional trace 
    if (return_trace) {
      p_stop_true_next <- if (k < K) {
        p_stop_next_given_theta(theta_true, z_stat[k], I_cum, I_cum + Ip_stage, a_crit[k+1], b_crit[k+1])
      } else NA_real_
      
      trace <- bind_rows(trace, tibble(
        k = k,
        theta_true = theta_true,
        r = r[k],
        I_cum = I_cum,
        z = z_stat[k],
        mu_post = mu_post[k],
        tau_post = tau_post[k],
        p_stop_true_next = p_stop_true_next,
        stopped_here = ifelse(!is.na(stop_k) && stop_k == k, TRUE, FALSE)
      ))
    }
    
    if (!is.na(stop_k)) break
    if (k < K) {
      r[k + 1] <- r_next_wp(mu_post[k], tau_post[k], a, delta, r_bounds)
    }
  }
  
  if (is.na(stop_k)) stop_k <- K
  
  if (loss_mode == "accum") {
    final_post <- sum(loss_post_stage[1:stop_k])
    final_real <- sum(loss_real_stage[1:stop_k])
  } else {
    I1_tot <- sum(I1_stage[1:stop_k])
    I2_tot <- sum(I2_stage[1:stop_k])
    final_post <- expected_loss_mu_tau_I(mu_post[stop_k], tau_post[stop_k], I1_tot, I2_tot, a, delta)
    final_real <- loss_function(theta_true, I1_tot, I2_tot, a, delta)
  }
  out <- list(
    final_post_loss = final_post,
    final_realised_loss = final_real,
    stop_k = stop_k
  )
  if (return_trace) out$trace <- trace
  out
}


run_grouped_correct_chunked <- function(
    cl, priors, nsims,
    loss_mode = c("accum","terminal"),
    chunk_size = 100000L,
    K, sigma, I_theta_fix, delta,
    r_bounds = c(0.05, 20)
) {
  loss_mode <- match.arg(loss_mode)
  Z <- nrow(priors)
  n_chunks <- ceiling(nsims / chunk_size)
  
  est_post <- est_real <- numeric(Z)
  mcse_post <- mcse_real <- numeric(Z)
  
  pb <- txtProgressBar(min = 0, max = Z * n_chunks, style = 3)
  step <- 0L
  
  for (z in 1:Z) {
    theta_0 <- priors$theta_0[z]
    theta_tau_0 <- priors$theta_tau0[z]
    
    S_post <- Q_post <- 0
    S_real <- Q_real <- 0
    n_tot <- 0L
    
    for (c in 1:n_chunks) {
      offset <- (c - 1L) * chunk_size
      n_this <- min(chunk_size, nsims - offset)
      
      losses_list <- parLapply(
        cl, 1:n_this,
        function(i, z, offset, K, sigma, I_theta_fix, delta, theta_0, theta_tau_0, loss_mode, r_bounds) {
          set.seed(seed_ij(i, z, offset))
          sim <- simulate_BRAR_grp_thetapr_final_loss_golden(
            K = K, sigma = sigma, I_theta_fix = I_theta_fix, delta = delta,
            theta_0 = theta_0, theta_tau_0 = theta_tau_0,
            loss_mode = loss_mode,
            r_bounds = r_bounds,
            return_trace = FALSE
          )
          c(final_post_loss = sim$final_post_loss,
            final_realised_loss = sim$final_realised_loss)
        },
        z = z, offset = offset,
        K = K, sigma = sigma, I_theta_fix = I_theta_fix, delta = delta,
        theta_0 = theta_0, theta_tau_0 = theta_tau_0,
        loss_mode = loss_mode, r_bounds = r_bounds
      )
      
      mat <- do.call(rbind, losses_list)
      post <- mat[, "final_post_loss"]
      real <- mat[, "final_realised_loss"]
      
      S_post <- S_post + sum(post); Q_post <- Q_post + sum(post * post)
      S_real <- S_real + sum(real); Q_real <- Q_real + sum(real * real)
      n_tot  <- n_tot + length(post)
      
      step <- step + 1L
      setTxtProgressBar(pb, step)
    }
    
    est_post[z] <- S_post / n_tot
    est_real[z] <- S_real / n_tot
    
    s2_post <- (Q_post - (S_post * S_post) / n_tot) / (n_tot - 1L)
    s2_real <- (Q_real - (S_real * S_real) / n_tot) / (n_tot - 1L)
    
    mcse_post[z] <- sqrt(s2_post / n_tot)
    mcse_real[z] <- sqrt(s2_real / n_tot)
  }
  
  close(pb)
  
  data.frame(
    prior = 1:nrow(priors),
    theta_0 = priors$theta_0,
    theta_tau0 = priors$theta_tau0,
    loss_mode = loss_mode,
    est_post = est_post, mcse_post = mcse_post,
    est_real = est_real, mcse_real = mcse_real
  )
}


num_cores <- max(1L, detectCores() - 1L)
cl <- makeCluster(num_cores)

clusterExport(
  cl,
  varlist = c(
    "seed_ij","BASE_SEED","STRIDE",
    "a","delta","inflation_factor","a_crit","b_crit",
    "loss_function","expected_loss_mu_tau_I",
    "r_next_wp","p_stop_next_given_theta",
    "simulate_BRAR_grp_thetapr_final_loss_golden","run_grouped_correct_chunked",
    "priors","nsims","I_theta_fix","K","sigma"
  ),
  envir = environment()
)

start.time <- Sys.time()
res <- run_grouped_correct_chunked(
  cl = cl, priors = priors, nsims = nsims,
  loss_mode = LOSS_MODE,
  chunk_size = 100000L,
  K = K, sigma = sigma,
  I_theta_fix = I_theta_fix, delta = delta,
  r_bounds = c(0.05, 20)
)
end.time <- Sys.time()

stopCluster(cl)
print(end.time - start.time)
print(res)


one <- simulate_BRAR_grp_thetapr_final_loss_golden(
  K = K, sigma = sigma, I_theta_fix = I_theta_fix, delta = delta,
  theta_0 = priors$theta_0[1], theta_tau_0 = priors$theta_tau0[1],
  loss_mode = LOSS_MODE,
  return_trace = TRUE
)
print(one$trace)