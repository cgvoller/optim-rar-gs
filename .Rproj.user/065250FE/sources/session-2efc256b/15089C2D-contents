# Create arrays using the following:
#         array_names: Vector of names to create arrays
#         nsims: Number of simulations
#         theta: Vector of thetas
#         Ks: How many analyses to store results for ()
initialise_arrays <- function(names2D = character(0),
                              names3D = character(0),
                              nsims,
                              theta,
                              Ks = NULL,
                              default_K = 5) {
  array_names <- c(names2D, names3D)
  
  if (length(names3D) > 0) {
    if (is.null(Ks)) {
      Ks_used <- rep(default_K, length(names3D))
    } else {
      stopifnot(length(names3D) == length(Ks))
      Ks_used <- Ks
    }
  } else {
    Ks_used <- numeric(0)
  }
  
  results <- setNames(lapply(seq_along(array_names), function(i) {
    name <- array_names[i]
    
    if (name %in% names2D) {
      array(0, dim = c(nsims, length(theta)))
    } else if (name %in% names3D) {
      K <- Ks_used[which(names3D == name)]
      array(0, dim = c(nsims, K, length(theta)))
    } else {
      array(0, dim = c(nsims, default_K, length(theta)))
    }
  }), array_names)
  
  return(results)
}

# Create arrays for storing results using bayesian methods
# Additional dimension to handle n priors
bayes_arrays <- function(nsims, theta, K,priors) {
  results <- list(
    n1_res = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    n2_res = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    #informationequal = array(0, dim = c(nsims, length(theta))),
    theta_hat = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    theta_hat_direct = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    tau1_res = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    tau2_res = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    mu1_res = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    mu2_res = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    x_1_res = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    x_2_res = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    ratio_res = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    prior_used = vector("list", nrow(priors))
  )
  return(results)
}

initialise_arrays_bayes <- function(names4D = character(0),
                                    names3D = character(0),
                                    names1D = character(0),
                                    nsims,
                                    theta,
                                    K,
                                    priors) {
  n_theta <- length(theta)
  n_priors <- nrow(priors)
  
  results4D <- lapply(names4D, \(x) array(0, dim = c(nsims, K, n_theta, n_priors)))
  names(results4D) <- names4D
  
  results3D <- lapply(names3D, \(x) array(0, dim = c(nsims, n_theta, n_priors)))
  names(results3D) <- names3D
  
  results1D <- lapply(names1D, \(x) vector("list", n_priors))
  names(results1D) <- names1D
  
  results <- c(results4D, results3D, results1D)
  return(results)
}


bayes_arrays_theta <- function(nsims, theta, K,priors) {
  results <- list(
    n1_res = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    n2_res = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    #informationequal = array(0, dim = c(nsims, length(theta))),
    theta_hat = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    theta_hat_tau = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    theta_hat_n = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    x_1_res = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    x_2_res = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    ratio_res = array(0, dim = c(nsims, K, length(theta),nrow(priors))),
    prior_used = vector("list", nrow(priors))
  )
  return(results)
}


