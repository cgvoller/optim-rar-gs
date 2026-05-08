## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/optim-rar-gs
##
## Purpose of script: Table 5: Fixed 3:2 sampling ratio
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
# No early stopping ------------------------------------------------------------
message("Simulate 3:2 Sampling - No early stopping")

simulate_fixed <- function(N, K, mu_1, mu_2, sigma, I_theta_fix,delta) {
  # Ratio, estimates of mu, theta
  muhat_1 <- muhat_2 <- theta_hat <- ratio <- numeric(K)
  # data from trial
  x_1 <- x_2 <- numeric(0)
  # Number of patients per arm
  n1 <- n2 <- n1.new <- n2.new <- numeric(K)
  for (k in 1:K) {
    if (k == 1) {
      # Initial equal allocation
      ratio[k] <- 3/2
      # Enroll patients
      n1.new[k] <- n2.new[k] <- N / K
      n1[k] <- n2[k] <- N / K
    } else {
      # Define ratio based on formula with previous estimate of theta
      ratio[k] <-  3/2
      # Patients in N1
      n1[k] <- (k / K) * (sigma ^ 2) * I_theta_fix * (1 + ratio[k])
      # Patient in N2
      n2[k] <- (k / K) * (sigma ^ 2) * I_theta_fix * (1 + (1 / ratio[k]))
      # Difference in patients from previous
      n1.new[k] <- n1[k] - n1[k - 1]
      n2.new[k] <- n2[k] - n2[k - 1]
      if(n1.new[k]<=0){
        n1[k]= n1[k - 1]
        n1.new[k] = n1[k] - n1[k - 1]
        #n2[k] = (((k/5)*I_theta_fix)^(-1)-1/(n2/sigma^2))^(-1)
        n2[k] = (k*n1[k]*sigma^2*I_theta_fix)/(5*n1[k]-k*sigma^2*I_theta_fix)
        n2.new[k] = n2[k]-n2[k-1]
      }
      if(n2.new[k]<=0){
        n2[k]= n2[k - 1]
        n2.new[k] = n2[k] - n2[k - 1]
        #n1[k] = (((k/5)*I_theta_fix)^(-1)-1/(n2[k]/sigma^2))^(-1)
        n1[k] = (k*n2[k]*sigma^2*I_theta_fix)/(5*n2[k]-k*sigma^2*I_theta_fix)
        n1.new[k] = n1[k]-n1[k-1]
      }
    }
    # Sample once from normal distribution
    x_1[k] <- rnorm(1, mu_1, sqrt(sigma^2/n1.new[k]))
    x_2[k] <- rnorm(1, mu_2, sqrt(sigma^2/n2.new[k]))
    
    muhat_1[k] <- sum(n1.new[1:k]*x_1[1:k])/sum(n1.new)
    muhat_2[k] <- sum(n2.new[1:k]*x_2[1:k])/sum(n2.new)
    # Treatment effect, difference of two means
    theta_hat[k] <- muhat_1[k] - muhat_2[k]
  }
  return(list(n1 = n1, n2 = n2, theta_hat = theta_hat))
}

# Initialise arrays
fptb5_results <- initialise_arrays(
  names3D =  c("theta_hats"),
  names2D = c("fptb5_EN1", "fptb5_EN2"),
  nsims = JenNsims,
  theta = theta,
  Ks=K
)

# ---------- Setup cluster ----------
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterSetRNGStream(cl, iseed = 42)
clusterExport(cl,
              varlist = c(
                "simulate_fixed",
                "I_theta_fix",
                "delta",
                "theta",
                "JenNsims",
                "a",
                "K"
              ))

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
  clusterExport(cl,varlist=c("j"))
  simresults <- parLapply(cl,1:JenNsims, \(x) {
    # After each run of theta, print as progress
    if (x %% JenNsims == 00)
      print(x)
    simulate_fixed(
      N = 100,
      K = K,
      mu_1 = theta[j],
      mu_2 = 0,
      sigma = 1,
      I_theta_fix = I_theta_fix,
      delta=delta
    )
  })
  update_progress(JenNsims)
  for (i in 1:nsims) {
    fptb5_results$fptb5_EN1[i, j] = unlist(simresults[[i]]$n1[K])
    fptb5_results$fptb5_EN2[i, j] = unlist(simresults[[i]]$n2[K])
    fptb5_results$theta_hats[i,,j] <- unlist(simresults[[i]]$theta_hat)
  }
}

stopCluster(cl)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# Calculate pooled expected sample sizes
EN1_pool_se =  apply(fptb5_results$fptb5_EN1, 2, \(x) se_funct2(x))
EN2_pool_se = apply(fptb5_results$fptb5_EN2, 2, \(x) se_funct2(x))
total_pool = round2(EN1_pool_se["mean", ] + EN2_pool_se["mean", ],digit=dig)
EN1_pool = apply(EN1_pool_se, 2, \(x) tidy_se(x, digit = dig))
EN2_pool = apply(EN2_pool_se, 2, \(x) tidy_se(x, digit = dig))


# Early stopping ---------------------------------------------------------------
message("Simulate 3:2 Sampling - With early stopping")

simulate_fixed_early <- function(N, K, mu_1, mu_2, sigma, I_theta_fix,delta) {
  # Estimates of mu, theta and ratio
  muhat_1 <- muhat_2 <- theta_hat <- ratio <- numeric(K)
  # data from trial
  x_1 <- x_2 <- numeric(0)
  # Information
  I1 <- I2 <- Ip <- numeric(K)
  # Total and new patients per arm at analysis K
  n1 <- n2 <- n1.new <- n2.new <- numeric(K)
  # Early stopping vars
  z <- numeric(K)
  # Inflation factor for GSD
  I_max <- I_theta_fix * inflation_factor
  for (k in 1:K) {
    if (k == 1) {
      # Initial equal allocation
      ratio[k] <- 3/2
      # Enroll patients
      n1.new[k] <- n1[k] <- (k / K) * (sigma ^ 2) * I_max * (1 + ratio[k])
      n2[k] <- n2.new[k] <- (k / K) * (sigma ^ 2) * I_max *(1 + (1 / ratio[k]))
      
    } else {
      # Define ratio based on formula with previous estimate of theta
      ratio[k] <- 3/2
      # Patients in N1
      n1[k] <- (k / K) * (sigma ^ 2) * I_max * (1 + ratio[k])
      # Patient in N2
      n2[k] <- (k / K) * (sigma ^ 2) * I_max *(1 + (1 / ratio[k]))
      # Difference in patients from previous
      n1.new[k] <- n1[k] - n1[k - 1]
      n2.new[k] <- n2[k] - n2[k - 1]
      if(n1.new[k]<=0){
        n1[k]= n1[k - 1]+0.0001
        n1.new[k] = 0
        n2[k] = (k*n1[k]*sigma^2*I_max)/(5*n1[k]-k*sigma^2*I_max)
        n2.new[k] = n2[k]-n2[k-1]
      }
      if(n2.new[k]<=0){
        n2[k]= n2[k - 1]+0.0001
        n2.new[k] = 0
        n1[k] = (k*n2[k]*sigma^2*I_max)/(5*n2[k]-k*sigma^2*I_max)
        n1.new[k] = n1[k]-n1[k-1]
      }
    }
    # Sample once from normal distribution
    x_1[k] <- rnorm(1, mu_1, sqrt(sigma^2/n1.new[k]))
    x_2[k] <- rnorm(1, mu_2, sqrt(sigma^2/n2.new[k]))
    
    muhat_1[k] <- sum(n1.new[1:k]*x_1[1:k])/sum(n1.new)
    muhat_2[k] <- sum(n2.new[1:k]*x_2[1:k])/sum(n2.new)
    # Treatment effect, difference of two means
    theta_hat[k] <- muhat_1[k] - muhat_2[k]
    I1[k] <- n1[k] / sigma^2
    I2[k] <- n2[k] / sigma^2
    Ip[k] <- 1 / (1 / I1[k] + 1 / I2[k])
    z[k] <- (theta_hat[k]) * sqrt(Ip[k])
    # Early stopping (futility and efficacy boundaries)
    if (z[k] < a_crit[k] | z[k] > b_crit[k]) {
      break
    }
  }
  return(list(n1 = n1, n2 = n2, theta_hat = theta_hat, z = z))
}

# Initialise arrays
festb5_results <- initialise_arrays(
  names3D =  c("theta_hats","z"),
  names2D = c("festb5_EN1", "festb5_EN2"),
  nsims = nsims,
  theta = theta,
  Ks=c(K,K)
)

# Parallel coding
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterSetRNGStream(cl, iseed = 42)
clusterExport(cl,
              varlist = c(
                "simulate_fixed_early",
                "I_theta_fix",
                "delta",
                "theta",
                "JenNsims",
                "a",
                "K",
                "inflation_factor",
                "a_crit",
                "b_crit"
              ))

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
  clusterExport(cl,varlist=c("j"))
  simresults <- parLapply(cl, 1:JenNsims, \(x) {
    # After each run of theta, print as progress
    if (x %% JenNsims == 00)
      print(x)
    simulate_fixed_early(
      N = 100,
      K = K,
      mu_1 = theta[j],
      mu_2 = 0,
      sigma = 1,
      I_theta_fix = I_theta_fix,
      delta=delta
    )
  })
  update_progress(JenNsims)
  for (i in 1:nsims) {
    festb5_results$festb5_EN1[i,j] = unlist(simresults[[i]]$n1[max(which(simresults[[i]]$n1!=0))])
    festb5_results$festb5_EN2[i,j] = unlist(simresults[[i]]$n2[max(which(simresults[[i]]$n2!=0))])
    festb5_results$theta_hats[i,,j] <- unlist(simresults[[i]]$theta_hat)
    festb5_results$z[i,,j] <- unlist(simresults[[i]]$z)
  }
}
stopCluster(cl)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


# Calculate pooled expected sample sizes
EN1_pooles =  apply(festb5_results$festb5_EN1, 2, \(x) se_funct2(x))
EN2_pooles = apply(festb5_results$festb5_EN2, 2, \(x) se_funct2(x))
total_pooles = round2(EN1_pooles["mean", ] + EN2_pooles["mean", ],digit=dig)
EN1_pooles = apply(EN1_pooles, 2, \(x) tidy_se(x, digit = dig))
EN2_pooles = apply(EN2_pooles, 2, \(x) tidy_se(x, digit = dig))

EN1_pooles
EN2_pooles

# Early stopping + Delay -------------------------------------------------------
message("Simulate 3:2 Sampling - With early stopping and delay")

simulate_fixed_earlydl <- function(N, K, mu_1, mu_2, sigma, I_theta_fix,delta) {
  # Estimates of mu, theta and ratio
  muhat_1 <- muhat_2 <- theta_hat <- ratio <- numeric(K)
  # data from trial
  x_1 <- x_2 <- numeric(0)
  # Information
  I1 <- I2 <- Ip <- numeric(K)
  # Total and new patients per arm at analysis K
  n1 <- n2 <- n1.new <- n2.new <- numeric(K)
  # Early stopping vars
  z <- numeric(K)
  # Inflation factor for GSD
  I_max <- I_theta_fix * inflation_factor
  for (k in 1:K) {
    if (k < 3) {
      # Initial equal allocation
      ratio[k] <- 3 / 2
      # Enroll patients
      n1.new[k] <- (1 / K) * (sigma^2) * I_max * (1 + ratio[k])
      n1[k] <- (k / K) * (sigma^2) * I_max * (1 + ratio[k])
      n2.new[k] <- (1 / K) * (sigma^2) * I_max * (1 + (1 / ratio[k]))
      n2[k] <- (k / K) * (sigma^2) * I_max * (1 + (1 / ratio[k]))
    } else {
      # Define ratio based on formula with previous estimate of theta
      ratio[k] <-  3/2
      # Patients in N1
      n1[k] <- (k / K) * (sigma ^ 2) * I_max * (1 + ratio[k])
      # Patient in N2
      n2[k] <- (k / K) * (sigma ^ 2) * I_max *(1 + (1 / ratio[k]))
      # Difference in patients from previous
      n1.new[k] <- n1[k] - n1[k - 1]
      n2.new[k] <- n2[k] - n2[k - 1]
      if(n1.new[k]<=0){
        n1[k]= n1[k - 1]
        n1.new[k] = 0
        n2[k] = (k*n1[k]*sigma^2*I_max)/(5*n1[k]-k*sigma^2*I_max)
        n2.new[k] = n2[k]-n2[k-1]
      }
      if(n2.new[k]<=0){
        n2[k]= n2[k - 1]
        n2.new[k] = 0
        n1[k] = (k*n2[k]*sigma^2*I_max)/(5*n2[k]-k*sigma^2*I_max)
        n1.new[k] = n1[k]-n1[k-1]
      }
    }
    # Sample once from normal distribution
    x_1[k] <- rnorm(1, mu_1, sqrt(sigma^2/n1.new[k]))
    x_2[k] <- rnorm(1, mu_2, sqrt(sigma^2/n2.new[k]))
    
    muhat_1[k] <- sum(n1.new[1:k]*x_1[1:k])/sum(n1.new)
    muhat_2[k] <- sum(n2.new[1:k]*x_2[1:k])/sum(n2.new)
    # Treatment effect, difference of two means
    theta_hat[k] <- muhat_1[k] - muhat_2[k]
    I1[k] <- n1[k] / sigma^2
    I2[k] <- n2[k] / sigma^2
    Ip[k] <- 1 / (1 / I1[k] + 1 / I2[k])
    z[k] <- (theta_hat[k]) * sqrt(Ip[k])
    if (k>1){
      if (z[k-1] < a_crit[k-1] | z[k-1] > b_crit[k-1]) {
        break
      }
    }
  }
  return(list(n1 = n1, n2 = n2, theta_hat = theta_hat, z = z))
}

# Initialise arrays
fesdtb5_results <- initialise_arrays(
  names3D =  c("theta_hats","z"),
  names2D = c("fesdtb5_EN1", "fesdtb5_EN2"),
  nsims = JenNsims,
  theta = theta,
  Ks=c(K,K)
)

# Parallel coding
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterSetRNGStream(cl, iseed = 42)
clusterExport(cl,
              varlist = c(
                "simulate_fixed_earlydl",
                "I_theta_fix",
                "delta",
                "theta",
                "JenNsims",
                "a",
                "K",
                "inflation_factor",
                "a_crit",
                "b_crit"
              ))

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
# For each theta, loop through RAR simulation nsims times
for (j in 1:length(theta)) {
  clusterExport(cl,varlist=c("j"))
  simresults <- parLapply(cl, 1:JenNsims, \(x) {
    # After each run of theta, print as progress
    if (x %% JenNsims == 00)
      print(x)
    simulate_fixed_earlydl(
      N = 100,
      K = 5,
      mu_1 = theta[j],
      mu_2 = 0,
      sigma = 1,
      I_theta_fix = I_theta_fix,
      delta=delta
    )
  })
  update_progress(JenNsims)
  for (i in 1:nsims) {
    fesdtb5_results$fesdtb5_EN1[i,j] = unlist(simresults[[i]]$n1[max(which(simresults[[i]]$n1!=0))])
    fesdtb5_results$fesdtb5_EN2[i,j] = unlist(simresults[[i]]$n2[max(which(simresults[[i]]$n2!=0))])
    fesdtb5_results$theta_hats[i,,j] <- unlist(simresults[[i]]$theta_hat)
    fesdtb5_results$z[i,,j] <- unlist(simresults[[i]]$z)
  }
}
stopCluster(cl)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


# Calculate pooled expected sample sizes
EN1_poolesd =  apply(fesdtb5_results$fesdtb5_EN1, 2, \(x) se_funct2(x))
EN2_poolesd = apply(fesdtb5_results$fesdtb5_EN2, 2, \(x) se_funct2(x))
total_poolesd = round2(EN1_poolesd["mean", ] + EN2_poolesd["mean", ],digit=dig)
EN1_poolesd = apply(EN1_poolesd, 2, \(x) tidy_se(x, digit = dig))
EN2_poolesd = apply(EN2_poolesd, 2, \(x) tidy_se(x, digit = dig))


# Output -----------------------------------------------------------------------

# Table
freq_tb5 <- create_tab_fun(
  EN1_pool,
  EN2_pool,
  total_pool,
  EN1_pooles,
  EN2_pooles,
  total_pooles,
  EN1_poolesd,
  EN2_poolesd,
  total_poolesd,
  headertxt = "Designs with fixed 3:2 sampling ratio",
  footnotetxt = paste("mean (standard error) to", dig, "decimal place")
)

freq_tb5
# Remove objects ---------------------------------------------------------------

rm(
  EN1_pool,
  EN2_pool,
  total_pool,
  EN1_pooles,
  EN2_pooles,
  total_pooles,
  EN1_poolesd,
  EN2_poolesd,
  total_poolesd
)