## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/BRAR_gsd
##
## Purpose of script: Calculate posterior distributions given norm-norm
##
## Author: Corey Voller
##
## Date Created: 14-03-2025
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
## Normal Posterior Mean -------------------------------------------------------

posterior_mean <- function(prior_mean, prior_tau, data, tau, n,k) {
  result <- tryCatch({
    # Check if the lengths of n and data are compatible
    if (length(n) != length(data)) {
      stop("Length of n and data must be the same. Found n of length ", length(n), " and data of length ", length(data))
    }
    # Calculate post mean
    posterior_mean_result <- ((prior_tau) * prior_mean + sum(n * (tau) * (data))) / (sum(n* tau) + prior_tau)
    return(posterior_mean_result)
    
  }, error = function(e) {
    # Catch the error and return a message
    cat("Error in posterior_mean (K=):",k, e$message, "\n")
    cat("n: ",paste(n),"\n")
    cat("data: ",paste(data),"\n")
    return(NA)  # Return NA or any default value in case of error
  })
  
  return(result)
}

## Normal Posterior tau --------------------------------------------------------

posterior_tau <- function(prior_tau,n,tau){
  prior_tau + sum(n*tau)
}

# Normal Posterior Theta -------------------------------------------------------

# posterior_theta <- function(prior_mean, prior_tau, data, tau, n1, n2) {
#   result <- tryCatch({
#     # Check if the lengths of n and data are compatible
#     if (length(n1) != length(data) || length(n2) != length(data)) {
#       stop(
#         paste0(
#           "Length mismatch: n1 = ", length(n1),
#           ", n2 = ", length(n2),
#           ", data = ", length(data)
#         )
#       )
#     }
#     # Assume equal precision
#     lh <- (n1 * tau * n2 * tau) / (n1*tau + n2*tau)
#     posterior_precision <- prior_tau + sum(lh)
#     # Calculate posterior mean for theta
#     posterior_mean_result <- (prior_tau * prior_mean + sum(data*lh))/(posterior_precision)
#     return(posterior_mean_result)
#   }, error = function(e) {
#     # Catch the error and return a message
#     cat("Error in posterior_theta:", e$message, "\n")
#     return(NA)  # Return NA or any default value in case of error
#   })
#   
#   return(result)
# }
# posterior_theta <- function(prior_mean, prior_tau, data, tau, n1, n2) {
#   result <- tryCatch({
#     # Check if the lengths of n and data are compatible
#     if (length(n1) != length(data) || length(n2) != length(data)) {
#       stop(
#         paste0(
#           "Length mismatch: n1 = ", length(n1),
#           ", n2 = ", length(n2),
#           ", data = ", length(data)
#         )
#       )
#     }
#     # Assume equal precision
#     posterior_precision <- prior_tau + sum((n1*tau*n2*tau)/(tau*n1+n2*tau))
#     # Calculate posterior mean for theta
#     posterior_mean_result <- (prior_tau * prior_mean + sum(data*((n1*tau*n2*tau)/(tau*n1+n2*tau))))/(posterior_precision)
#     return(posterior_mean_result)
#   }, error = function(e) {
#     # Catch the error and return a message
#     cat("Error in posterior_theta:", e$message, "\n")
#     return(NA)  # Return NA or any default value in case of error
#   })
#   
#   return(result)
# }


posterior_theta <- function(prior_mean, prior_tau, data, tau, n1, n2) {
  result <- tryCatch({
    # Check if the lengths of n and data are compatible
    if (length(n1) != length(data) || length(n2) != length(data)) {
      stop(
        paste0(
          "Length mismatch: n1 = ", length(n1),
          ", n2 = ", length(n2),
          ", data = ", length(data)
        )
      )
    }
    # Assume equal precision
    posterior_precision <- prior_tau +sum(prior_tau + sum(1/(2/(tau))))
    # Calculate posterior mean for theta
    posterior_mean_result <- (prior_tau * prior_mean + sum(data*(tau*(1/n1 + 1/n2))))/(posterior_precision)
    return(posterior_mean_result)
  }, error = function(e) {
    # Catch the error and return a message
    cat("Error in posterior_theta:", e$message, "\n")
    return(NA)  # Return NA or any default value in case of error
  })
  
  return(result)
}

