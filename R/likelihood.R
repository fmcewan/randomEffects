# Validation method 
validate_method <- function(method) {
  
  if (!is.character(method) || is.na(method)) {
    stop("Invalid method name.")  
  }

}

# cluster_loglik and marginal_loglik
cluster_loglik <- function(theta, cd, method, K, draws) {

  methods <- c("gh", "laplace", "agh", "is")

  validate_method(method)

  result <- switch(
    method, 
    "gh" = gauss_hermite(theta, cd, K),
    "laplace" = laplace(theta, cd),
    "agh" = adaptive_gauss_hermite(theta, cd, K),
    "is" = importance_sampling(theta, cd, draws),
    stop("Method not yet implemented.")
  )

  return(result)

}

