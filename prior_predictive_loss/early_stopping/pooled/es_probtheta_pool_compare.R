############################################################
## POOLED: Compare Method (iii) vs Method (iv)
##
##  r3 = whole posterior closed-form ratio (no stop weighting)
##  r4 = golden-section search with stop probability weighting
##
## Runs in the same chunked parallel structure as your scripts.
## Also samples r3/r4 ratios and saves them for plotting/tables.
############################################################

if (RUN_POOLED) {
  
  message("POOLED: Compare r3 (whole posterior closed form) vs r4 (stop-weighted golden search)")
  
  # ----------------------------
  # 1) Loss functions
  # ----------------------------
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
  
  loss_theta_I <- function(theta, I1, I2, w_exp) {
    ifelse(theta > 0,
           I1 + I2 * exp(w_exp * theta),
           I2 + I1 * exp(-w_exp * theta))
  }
  
  # ----------------------------
  # 2) Method (iii): whole posterior closed-form ratio r3
  # ----------------------------
  ratio_whole_posterior_closed_form <- function(mu_post, tau_post, a, delta,
                                                r_bounds = c(0.05, 20)) {
    w <- log(a) / delta
    v <- 1 / tau_post
    sd <- sqrt(v)
    
    p_pos <- pnorm(mu_post / sd)
    p_neg <- 1 - p_pos
    
    C2_pos <- exp(w * mu_post + 0.5 * w^2 * v) * pnorm((mu_post + w * v) / sd)
    C1_neg <- exp(-w * mu_post + 0.5 * w^2 * v) * pnorm((w * v - mu_post) / sd)
    
    C1 <- p_pos + C1_neg
    C2 <- C2_pos + p_neg
    
    r <- sqrt(C2 / C1)
    max(r_bounds[1], min(r_bounds[2], r))
  }
  
  # ----------------------------
  # 3) One-step stop probability (Normal tails)
  #     IMPORTANT: I_prev and I_k are cumulative pooled infos (I_{θ,k-1}, I_{θ,k})
  # ----------------------------
  p_stop_at_k_given_theta <- function(theta, z_prev, I_prev, I_k, a_k, b_k) {
    dI <- I_k - I_prev
    if (dI <= 0) return(0)
    
    m <- (z_prev * sqrt(I_prev) + theta * dI) / sqrt(I_k)
    s <- sqrt(dI / I_k)
    
    pnorm((a_k - m) / s) + (1 - pnorm((b_k - m) / s))
  }
  
  stop_weighted_integral_mc <- function(theta_draws,
                                        z_prev, I_prev, I_k,
                                        a_k, b_k,
                                        I1_k, I2_k,
                                        w_exp) {
    pS <- p_stop_at_k_given_theta(theta_draws, z_prev, I_prev, I_k, a_k, b_k)
    L  <- loss_theta_I(theta_draws, I1_k, I2_k, w_exp)
    mean(pS * L)
  }
  
  # ----------------------------
  # 4) Method (iv): golden-section / optimize ratio r4
  # ----------------------------
  choose_ratio_golden <- function(mu_prev, tau_prev,
                                  z_prev, I_prev,
                                  k, K, I_max, sigma,
                                  a_crit, b_crit,
                                  a, delta,
                                  M_theta = 10000L,
                                  r_bounds = c(0.05, 20),
                                  seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    
    theta_draws <- rnorm(M_theta, mean = mu_prev, sd = 1 / sqrt(tau_prev))
    
    w_exp <- log(a) / delta
    a_k <- a_crit[k]
    b_k <- b_crit[k]
    
    f_logr <- function(logr) {
      r <- exp(logr)
      
      # Candidate cumulative n1,n2 at look k (pooled schedule)
      n1_k <- (k / K) * (sigma^2) * I_max * (1 + r)
      n2_k <- (k / K) * (sigma^2) * I_max * (1 + 1 / r)
      
      I1_k <- n1_k / sigma^2
      I2_k <- n2_k / sigma^2
      I_k  <- 1 / (1 / I1_k + 1 / I2_k)  # pooled info at look k
      
      stop_weighted_integral_mc(theta_draws,
                                z_prev = z_prev, I_prev = I_prev, I_k = I_k,
                                a_k = a_k, b_k = b_k,
                                I1_k = I1_k, I2_k = I2_k,
                                w_exp = w_exp)
    }
    
    opt <- optimize(f_logr, interval = log(r_bounds))
    exp(opt$minimum)
  }
  
  # ----------------------------
  # 5) One pooled replicate producing BOTH r3 and r4 sequences
  #    policy determines which ratio is actually used for allocation.
  # ----------------------------
  simulate_pooled_r3_r4 <- function(N, K, sigma, I_theta_fix, delta,
                                    theta_0, theta_tau_0,
                                    policy = c("r3", "r4"),
                                    M_theta = 10000L,
                                    r_bounds = c(0.05, 20)) {
    
    policy <- match.arg(policy)
    
    w_exp <- log(a) / delta
    I_max <- I_theta_fix * inflation_factor
    
    theta_true <- rnorm(1, mean = theta_0, sd = sqrt(1 / theta_tau_0))
    mu_1 <- theta_true
    mu_2 <- 0
    
    # store both ratio sequences
    r3 <- r4 <- ratio_used <- numeric(K)
    
    x1 <- x2 <- muhat1 <- muhat2 <- theta_hat <- numeric(K)
    I1 <- I2 <- Ip <- numeric(K)
    n1 <- n2 <- n1_new <- n2_new <- numeric(K)
    
    mu_post <- tau_post <- numeric(K)
    z_stat <- numeric(K)
    
    loss_post <- loss_real <- numeric(K)
    stop_k <- NA_integer_
    
    for (k in 1:K) {
      if (k == 1) {
        r3[k] <- 1
        r4[k] <- 1
        ratio_used[k] <- 1
        
        n1[k] <- n2[k] <- (k / K) * 2 * I_max
        n1_new[k] <- n1[k]
        n2_new[k] <- n2[k]
        
      } else {
        # compute BOTH candidate ratios from the same posterior state (k-1)
        r3[k] <- ratio_whole_posterior_closed_form(mu_post[k-1], tau_post[k-1], a, delta, r_bounds)
        
        r4[k] <- choose_ratio_golden(mu_post[k-1], tau_post[k-1],
                                     z_stat[k-1], Ip[k-1],
                                     k, K, I_max, sigma,
                                     a_crit, b_crit,
                                     a, delta,
                                     M_theta = M_theta,
                                     r_bounds = r_bounds)
        
        ratio_used[k] <- if (policy == "r3") r3[k] else r4[k]
        
        # cumulative n at look k implied by chosen ratio
        n1[k] <- (k / K) * (sigma^2) * I_max * (1 + ratio_used[k])
        n2[k] <- (k / K) * (sigma^2) * I_max * (1 + 1 / ratio_used[k])
        
        n1_new[k] <- n1[k] - n1[k-1]
        n2_new[k] <- n2[k] - n2[k-1]
        
        # guards to ensure positive increments
        if (n1_new[k] <= 0) {
          n1[k] <- n1[k-1] + 1e-5
          n1_new[k] <- n1[k] - n1[k-1]
          n2[k] <- (k * n1[k] * sigma^2 * I_max) / (K * n1[k] - k * sigma^2 * I_max)
          n2_new[k] <- n2[k] - n2[k-1]
        }
        if (n2_new[k] <= 0) {
          n2[k] <- n2[k-1] + 1e-5
          n2_new[k] <- n2[k] - n2[k-1]
          n1[k] <- (k * n2[k] * sigma^2 * I_max) / (K * n2[k] - k * sigma^2 * I_max)
          n1_new[k] <- n1[k] - n1[k-1]
        }
      }
      
      # stage means for new patients
      x1[k] <- rnorm(1, mean = mu_1, sd = sqrt(sigma^2 / n1_new[k]))
      x2[k] <- rnorm(1, mean = mu_2, sd = sqrt(sigma^2 / n2_new[k]))
      
      # cumulative means + theta_hat
      muhat1[k] <- sum(n1_new[1:k] * x1[1:k]) / sum(n1_new[1:k])
      muhat2[k] <- sum(n2_new[1:k] * x2[1:k]) / sum(n2_new[1:k])
      theta_hat[k] <- muhat1[k] - muhat2[k]
      
      # infos
      I1[k] <- n1[k] / sigma^2
      I2[k] <- n2[k] / sigma^2
      Ip[k] <- 1 / (1 / I1[k] + 1 / I2[k])
      
      # posterior update
      tau_post[k] <- theta_tau_0 + Ip[k]
      mu_post[k]  <- (theta_tau_0 * theta_0 + Ip[k] * theta_hat[k]) / tau_post[k]
      
      # losses at look k
      loss_post[k] <- expected_loss_mu_tau_I(mu_post[k], tau_post[k], I1[k], I2[k], a, delta)
      loss_real[k] <- loss_theta_I(theta_true, I1[k], I2[k], w_exp)
      
      # fixed stopping rule
      z_stat[k] <- theta_hat[k] * sqrt(Ip[k])
      if (z_stat[k] < a_crit[k] || z_stat[k] > b_crit[k]) { stop_k <- k; break }
    }
    
    if (is.na(stop_k)) stop_k <- K
    
    # return ratios + final losses
    c(setNames(r3, paste0("r3_", 1:K)),
      setNames(r4, paste0("r4_", 1:K)),
      final_post_loss = loss_post[stop_k],
      final_realised_loss = loss_real[stop_k])
  }
  
  # ----------------------------
  # 6) Chunked runner (parallel) over priors, like your other scripts
  #    Also saves sampled ratios to disk for plots/tables.
  # ----------------------------
  run_pooled_r3_r4_chunked <- function(
    cl,
    priors,
    nsims,
    policy = c("r3", "r4"),
    chunk_size = 100000L,
    N = 100L,
    K = 5L,
    sigma = 1,
    I_theta_fix,
    delta,
    M_theta = 10000L,
    r_bounds = c(0.05, 20),
    ratio_sample_frac = 0.01,
    ratio_out_dir = file.path("output", "ratio_samples_pooled")
  ) {
    policy <- match.arg(policy)
    dir.create(ratio_out_dir, showWarnings = FALSE, recursive = TRUE)
    
    Z <- nrow(priors)
    n_chunks <- ceiling(nsims / chunk_size)
    
    est_post <- est_real <- numeric(Z)
    mcse_post <- mcse_real <- numeric(Z)
    
    ratio_cols <- c(paste0("r3_", 1:K), paste0("r4_", 1:K))
    
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
          function(i, z, offset, N, K, sigma, I_theta_fix, delta,
                   theta_0, theta_tau_0, policy, M_theta, r_bounds) {
            set.seed(seed_ij(i, z, offset))
            simulate_pooled_r3_r4(
              N = N, K = K,
              sigma = sigma,
              I_theta_fix = I_theta_fix,
              delta = delta,
              theta_0 = theta_0,
              theta_tau_0 = theta_tau_0,
              policy = policy,
              M_theta = M_theta,
              r_bounds = r_bounds
            )
          },
          z = z, offset = offset,
          N = N, K = K, sigma = sigma,
          I_theta_fix = I_theta_fix, delta = delta,
          theta_0 = theta_0, theta_tau_0 = theta_tau_0,
          policy = policy, M_theta = M_theta, r_bounds = r_bounds
        )
        
        losses_mat <- do.call(rbind, losses_list)
        
        post <- losses_mat[, "final_post_loss"]
        real <- losses_mat[, "final_realised_loss"]
        
        S_post <- S_post + sum(post)
        Q_post <- Q_post + sum(post * post)
        
        S_real <- S_real + sum(real)
        Q_real <- Q_real + sum(real * real)
        
        n_tot <- n_tot + length(post)
        
        # --- save a sample of r3/r4 ratios for plotting ---
        if (ratio_sample_frac > 0) {
          keep_n <- max(1L, as.integer(ceiling(ratio_sample_frac * n_this)))
          keep_idx <- if (ratio_sample_frac >= 1) seq_len(n_this) else sample.int(n_this, keep_n)
          
          ratio_dt <- data.table::as.data.table(losses_mat[keep_idx, ratio_cols, drop = FALSE])
          ratio_dt[, prior := z]
          ratio_dt[, sim_id := offset + keep_idx]
          ratio_dt[, policy := policy]
          
          saveRDS(
            ratio_dt,
            file = file.path(ratio_out_dir, sprintf("policy_%s_prior_%03d_chunk_%05d.rds", policy, z, c))
          )
        }
        
        step <- step + 1L
        setTxtProgressBar(pb, step)
      }
      
      est_post[z] <- S_post / n_tot
      est_real[z] <- S_real / n_tot
      
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
  
  # ----------------------------
  # 7) Helpers: read ratio samples -> summary table + plot
  # ----------------------------
  summarise_ratio_samples <- function(ratio_out_dir, K) {
    files <- list.files(ratio_out_dir, pattern = "\\.rds$", full.names = TRUE)
    if (length(files) == 0) stop("No ratio sample files found in: ", ratio_out_dir)
    
    dt_list <- lapply(files, readRDS)
    dt <- data.table::rbindlist(dt_list, fill = TRUE)
    
    # long format
    ratio_cols <- c(paste0("r3_", 1:K), paste0("r4_", 1:K))
    dt_long <- data.table::melt(
      dt,
      id.vars = c("prior", "sim_id", "policy"),
      measure.vars = ratio_cols,
      variable.name = "method_stage",
      value.name = "ratio"
    )
    dt_long[, method := sub("_\\d+$", "", method_stage)]
    dt_long[, stage  := as.integer(sub("^.*_", "", method_stage))]
    
    # summary by prior/method/stage
    summ <- dt_long[, .(
      mean = mean(ratio, na.rm = TRUE),
      p10  = as.numeric(stats::quantile(ratio, 0.10, na.rm = TRUE)),
      p50  = as.numeric(stats::quantile(ratio, 0.50, na.rm = TRUE)),
      p90  = as.numeric(stats::quantile(ratio, 0.90, na.rm = TRUE))
    ), by = .(prior, method, stage)]
    
    list(dt_long = dt_long, summary = summ)
  }
  
  plot_ratio_boxplots <- function(dt_long) {
    # base R boxplot: method-stage
    dt_long$label <- paste0(dt_long$method, "_", dt_long$stage)
    boxplot(ratio ~ label, data = dt_long,
            las = 2, main = "Allocation ratios: r3 vs r4 (sampled)",
            ylab = "ratio")
    abline(h = 1, lty = 2)
  }
  
  # ----------------------------
  # 8) RUN EXAMPLE (same structure)
  # ----------------------------
  num_cores <- max(1L, parallel::detectCores() - 1L)
  cl <- parallel::makeCluster(num_cores)
  
  parallel::clusterExport(
    cl,
    varlist = c(
      "seed_ij", "BASE_SEED", "STRIDE",
      "expected_loss_mu_tau_I",
      "loss_theta_I",
      "ratio_whole_posterior_closed_form",
      "p_stop_at_k_given_theta",
      "stop_weighted_integral_mc",
      "choose_ratio_golden",
      "simulate_pooled_r3_r4",
      "run_pooled_r3_r4_chunked",
      "priors",
      "inflation_factor",
      "I_theta_fix", "delta", "a", "a_crit", "b_crit"
    ),
    envir = environment()
  )
  
  # run BOTH policies so you get comparable losses under r3 vs r4 allocation
  ratio_out_dir <- file.path("output", "ratio_samples_pooled")
  
  start.time <- Sys.time()
  res_r3 <- run_pooled_r3_r4_chunked(
    cl = cl,
    priors = priors,
    nsims = nsims,
    policy = "r3",
    chunk_size = chunk_size,
    N = 100L, K = K, sigma = sigma,
    I_theta_fix = I_theta_fix,
    delta = delta,
    M_theta = 10000L,
    r_bounds = c(0.05, 20),
    ratio_sample_frac = 0.01,
    ratio_out_dir = ratio_out_dir
  )
  
  res_r4 <- run_pooled_r3_r4_chunked(
    cl = cl,
    priors = priors,
    nsims = nsims,
    policy = "r4",
    chunk_size = chunk_size,
    N = 100L, K = K, sigma = sigma,
    I_theta_fix = I_theta_fix,
    delta = delta,
    M_theta = 10000L,
    r_bounds = c(0.05, 20),
    ratio_sample_frac = 0.01,
    ratio_out_dir = ratio_out_dir
  )
  
  end.time <- Sys.time()
  parallel::stopCluster(cl)
  print(end.time - start.time)
  
  # loss tables
  out_r3 <- data.frame(
    prior = 1:nrow(priors),
    theta_0 = priors$theta_0,
    theta_tau0 = priors$theta_tau0,
    policy = "r3",
    est_post = res_r3$est_post,
    mcse_post = res_r3$mcse_post,
    est_real = res_r3$est_real,
    mcse_real = res_r3$mcse_real
  )
  
  out_r4 <- data.frame(
    prior = 1:nrow(priors),
    theta_0 = priors$theta_0,
    theta_tau0 = priors$theta_tau0,
    policy = "r4",
    est_post = res_r4$est_post,
    mcse_post = res_r4$mcse_post,
    est_real = res_r4$est_real,
    mcse_real = res_r4$mcse_real
  )
  
  print(out_r3)
  print(out_r4)
  
  saveRDS(out_r3, file = file.path("output/expected_losses", "pooled_r3_loss.rds"))
  saveRDS(out_r4, file = file.path("output/expected_losses", "pooled_r4_loss.rds"))
  
  # ratio table + plot (from sampled files)
  ratio_summ <- summarise_ratio_samples(ratio_out_dir, K = K)
  print(ratio_summ$summary)
  plot_ratio_boxplots(ratio_summ$dt_long)
}
