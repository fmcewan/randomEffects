# Validation method 
validate_method <- function(method) {
  
  if (!is.character(method) || length(method) != 1 || is.na(method)) {
    stop("Invalid method name.")  
  }

}

# cluster_loglik
cluster_loglik <- function(theta, cd, method, K=NULL, draws=NULL) {

  validate_method(method)

  # TODO: Replace placeholders with the actual functions to call
  result <- switch(
    method, 
    "gh" = gauss_hermite(n),
    "laplace" = laplace(theta, cd),
    "agh" = adaptive_gauss_hermite(theta, cd, K),
    "is" = importance_sampling(theta, cd, draws),
    stop("Method not yet implemented.")
  )

  return(result)

}

# marginal_loglik
marginal_loglik <- function(theta, cd, method, K=NULL, draws=NULL) {
    sum(cluster_loglik(theta, cd, method, K, draws))
}

# fit_rint
