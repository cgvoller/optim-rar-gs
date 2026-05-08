## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/optim-rar-gs
##
## Purpose of script: Table 1: No early stopping (freq pooled+group)
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
message("Simulate RAR Pooled - No early stopping")
simulate_RAR <- function(K, mu_1, mu_2, sigma, I_theta_fix,delta) {
  # Ratio, estimates of mu, theta
  muhat_1 <- muhat_2 <- ratio <- numeric(K)
  theta_hat <- rep(NA, K)
  # data from trial
  x_1 <- x_2 <- numeric(0)
  # Number of patients per arm
  n1 <- n2 <- n1.new <- n2.new <- numeric(K)
  for (k in 1:K) {
    if (k == 1) {
      # Initial equal allocation
      ratio[k] <- 1
      # Enroll 20 patients when N = 100, K=5
      n1.new[k] <- n2.new[k] <- (k / K) * 2 * I_theta_fix
      n1[k] <- n2[k] <- (k / K) * 2 * I_theta_fix
    } else {
      # Define ratio based on formula with previous estimate of theta
      ratio[k] <-  a ** (theta_hat[k - 1]/ (2 * delta))
      # Patients in N1
      n1[k] <- (k / K) * (sigma ^ 2) * I_theta_fix * (1 + ratio[k])
      # Patients in N2
      n2[k] <- (k / K) * (sigma ^ 2) * I_theta_fix * (1 + (1 / ratio[k]))
      # Difference in patients from previous
      n1.new[k] <- n1[k] - n1[k - 1]
      n2.new[k] <- n2[k] - n2[k - 1]
      # Extreme cases lead to negative sample sizes, set n1 to small value and 
      # recalculate n2
      if(n1.new[k]<=0){
        n1[k]= n1[k - 1]+0.00001
        n1.new[k] = n1[k] - n1[k - 1]
        n2[k] = (k*n1[k]*sigma^2*I_theta_fix)/(5*n1[k]-k*sigma^2*I_theta_fix)
        n2.new[k] = n2[k]-n2[k-1]
      }
      if(n2.new[k]<=0){
        n2[k]= n2[k - 1]+0.00001
        n2.new[k] = n2[k] - n2[k - 1]
        n1[k] = (k*n2[k]*sigma^2*I_theta_fix)/(5*n2[k]-k*sigma^2*I_theta_fix)
        n1.new[k] = n1[k]-n1[k-1]
      }
    }
    # Sample once from normal distribution
    x_1[k] <- rnorm(1, mu_1, sqrt(sigma^2/n1.new[k]))
    x_2[k] <- rnorm(1, mu_2, sqrt(sigma^2/n2.new[k]))
    # Calculate treatment effect estimate
    muhat_1[k] <- sum(n1.new[1:k] * x_1[1:k]) / sum(n1.new)
    muhat_2[k] <- sum(n2.new[1:k] * x_2[1:k]) / sum(n2.new)
    theta_hat[k] <- muhat_1[k] - muhat_2[k]
  }
  return(list(n1 = n1, 
              n2 = n2, 
              theta_hat = theta_hat
              ))
}

# Initialise arrays
fptb1_results <- initialise_arrays(
  names3D =  c("theta_hats"),
  names2D = c("fptb1_EN1", "fptb1_EN2"),
  nsims = JenNsims,
  theta = theta,
  Ks=K
)

# ---------- Setup cluster ----------
num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)
clusterSetRNGStream(cl, iseed = 42)
clusterExport(
  cl, c(
    "simulate_RAR",
    "I_theta_fix",
    "delta",
    "theta",
    "JenNsims",
    "a",
    "K"
  ),
  envir=environment()
)

# ---------- Progress bar setup ----------
total_steps <- length(theta) * JenNsims
pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
progress_count <- 0
update_progress <- function(n = 1) {
  progress_count <<- progress_count + n
  setTxtProgressBar(pb, progress_count)
}


# Start timer
start.time <- Sys.time()
# For each theta, loop through RAR simulation JenNsims times
for (j in 1:length(theta)) {
  simresults <- parLapply(cl,1:JenNsims, \(x,theta_val) {
    simulate_RAR(
      K = K,
      mu_1 = theta_val,
      mu_2 = 0,
      sigma = 1,
      I_theta_fix = I_theta_fix,
      delta=delta
    )
  },theta_val=theta[j])
  update_progress(JenNsims)
  for (i in 1:JenNsims) {
    x <- simresults[[i]]
    fptb1_results$fptb1_EN1[i, j] = unlist(x$n1[K])
    fptb1_results$fptb1_EN2[i, j] = unlist(x$n2[K])
    fptb1_results$theta_hats[i,,j] = unlist(x$theta_hat)
  }
}
close(pb)
stopCluster(cl)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# Calculate pooled expected sample sizes
EN1_pool_se =  apply(fptb1_results$fptb1_EN1, 2, \(x) se_funct2(x))
EN2_pool_se = apply(fptb1_results$fptb1_EN2, 2, \(x) se_funct2(x))
total_pool = round2(EN1_pool_se["mean", ] + EN2_pool_se["mean", ],digit=dig)
EN1_pool = apply(EN1_pool_se, 2, \(x) tidy_se(x, digit = dig))
EN2_pool = apply(EN2_pool_se, 2, \(x) tidy_se(x, digit = dig))

# Simulate RAR Grouped ---------------------------------------------------------
message("Simulate RAR Grouped - No early stopping")

simulate_RAR_group <- function(K, mu_1, mu_2, sigma, I_theta_fix,delta) {
  muhat_1 <- muhat_2 <- theta_hat <- theta_hat_group <- ratio <- numeric(K)
  n1 <- n2 <- I1 <- I2 <- numeric(K)
  for (k in 1:K) {
    if (k == 1) {
      ratio[k] <- 1
      # Initial equal allocation
      n1[k] <- n2[k] <- (k / K) * (sigma ^ 2) * I_theta_fix * (1 + ratio[k])
    } else {
      ratio[k] <-  a ** (theta_hat_group[k - 1] / (2 * delta))
      n1[k] <- (1 / K) * (sigma ^ 2) * I_theta_fix * (1 + ratio[k])
      n2[k] <- (1 / K) * (sigma ^ 2) * I_theta_fix * (1 + (1 / ratio[k]))
    }
    # Sample once from normal distribution
    muhat_1[k] <- c(rnorm(1, mean = mu_1, sd = sqrt(sigma ^ 2 / n1[k])))
    muhat_2[k] <- c(rnorm(1, mean = mu_2, sd = sqrt((sigma ^ 2 / n2[k]))))
    # Calculate treatment effect estimate
    theta_hat[k] <- muhat_1[k] - muhat_2[k]
    theta_hat_group[k] = (1 / k) * sum(theta_hat[1:k])
    I1[k] <- n1[k] / sigma^2
    I2[k] <- n2[k] / sigma^2
  }
  return(list(
    n1 = n1,
    n2 = n2,
    theta_hat_group = theta_hat_group
  ))
}


# Initialise arrays
fgtb1_results <- initialise_arrays(
  names3D =  c("theta_group"),
  names2D = c("fgtb1_EN1", "fgtb1_EN2","final_loss"),
  nsims = JenNsims,
  theta = theta,
  Ks=K
)

# ---------- Setup cluster ----------
num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)
clusterSetRNGStream(cl, iseed = 42)
clusterExport(cl,
                c(
                "simulate_RAR_group",
                "loss_function",
                "I_theta_fix",
                "delta",
                "theta",
                "JenNsims",
                "a",
                "K"
              ),
              envir = environment()
              )

# ---------- Progress bar setup ----------
total_steps <- length(theta) * JenNsims
pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
progress_count <- 0
update_progress <- function(n = 1) {
  progress_count <<- progress_count + n
  setTxtProgressBar(pb, progress_count)
}

# Start timer
start.time <- Sys.time()

for (j in 1:length(theta)) {
  clusterExport(cl, varlist = c("j"))
  simresults <- parLapply(cl, 1:JenNsims, \(x,theta_val) {
    simulate_RAR_group(
      K = K,
      mu_1 = theta[j],
      mu_2 = 0,
      sigma = 1,
      I_theta_fix = I_theta_fix,
      delta = delta
    )
  },theta_val=theta[j])
  update_progress(JenNsims)
  for (i in 1:JenNsims) {
    fgtb1_results$fgtb1_EN1[i, j] = sum(simresults[[i]]$n1)
    fgtb1_results$fgtb1_EN2[i, j] = sum(simresults[[i]]$n2)
    fgtb1_results$theta_group[i, , j] = simresults[[i]]$theta_hat_group
  }
}
close(pb)
stopCluster(cl)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

EN1_group_se =  apply(fgtb1_results$fgtb1_EN1, 2, \(x) se_funct2(x))
EN2_group_se = apply(fgtb1_results$fgtb1_EN2, 2, \(x) se_funct2(x))
total_group = round2(EN1_group_se["mean", ] + EN2_group_se["mean", ],digit=dig)
EN1_group = apply(EN1_group_se, 2, \(x) tidy_se(x, digit = dig))
EN2_group = apply(EN2_group_se, 2, \(x) tidy_se(x, digit = dig))

# Target Values ----------------------------------------------------------------
EN1_target = round2(((1 + (a ** (theta / (
  2 * delta
)))) / 2) * 100 * sigma ^ 2,dig)
# Expected sample size for N2
EN2_target = round2(((1 + (a ** (-theta / (
  2 * delta
)))) / 2) * 100 * sigma ^ 2,dig)
target_total = round2(((1 + (a ** (theta / (
  2 * delta
)))) / 2) * 100 * sigma ^ 2 + ((1 + (a ** (-theta / (
  2 * delta
)))) / 2) * 100 * sigma ^ 2,dig)
# Output -----------------------------------------------------------------------

# Table
freq_tb1 <- create_tab_fun(EN1_target,
                         EN2_target,
                         target_total,
                         EN1_pool,
                         EN2_pool,
                         total_pool,
                         EN1_group,
                         EN2_group,
                         total_group,
                         headertxt = "Designs with no early stopping",
                         footnotetxt = paste("mean (standard error) to", dig, "decimal place"))
freq_tb1
# Remove objects ---------------------------------------------------------------

rm(
  EN1_pool,
  EN2_pool,
  total_pool,
  EN1_pool_se,
  EN2_pool_se,
  EN1_group,
  EN2_group,
  total_group,
  EN1_group_se,
  EN2_group_se
)
