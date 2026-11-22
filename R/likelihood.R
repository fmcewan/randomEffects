# Validation method 
validate_method <- function(method) {

  if (!is.character(method) || length(method) != 1 || is.na(method)) {
    stop("Invalid method name.")  
  }

}

#' Log-likelihood calculation for each cluster
#'
#' This method integrates out the random intercept of every cluster using 
#' the chosen method.
#'
#' @param theta Parameter vector (beta, log sigma).
#' @param cd Data prepared by cluster_data.
#' @param method Integration method: "gh", "laplace", "agh" or "is".
#' @param K Number of quadrature nodes (used by "gh" and "agh").
#' @param draws Importance-sampling draws from is_draws (used by "is").
#'
#' @return A numeric vector with one log-likelihood per cluster. For "is",
#'   it has an attribute mc_se, giving each cluster's Monte Carlo
#'   standard error.
#'
#' @export
#' @examples
#' data <- rint_sim(J = 50, m = 5, beta = c(-0.5, 1), sigma = 1)
#' cd <- cluster_data(data$y, data$X, factor(data$group))
#' cluster_loglik(c(-0.5, 1, 0), cd, method = "laplace")
cluster_loglik <- function(theta, cd, method, K=NULL, draws=NULL) {

  validate_method(method)

  switch(
    method, 
    "gh" = gh(theta, cd, K),
    "laplace" = laplace_loglik(theta, cd),
    "agh" = agh_loglik(theta, cd, K),
    "is" = importance_sampling(theta, cd, draws),
    stop("Method not yet implemented.")
  )

}

#' Marginal log-likelihood
#'
#' This method sums the cluster log-likelihoods from \code{\link{cluster_loglik}}.
#'
#' @inheritParams cluster_loglik
#'
#' @return The marginal log-likelihood (a single number).
#'
#' @export
#' @examples
#' data <- rint_sim(J = 50, m = 5, beta = c(-0.5, 1), sigma = 1)
#' cd <- cluster_data(data$y, data$X, factor(data$group))
#' marginal_loglik(c(-0.5, 1, 0), cd, method = "laplace")
marginal_loglik <- function(theta, cd, method, K=NULL, draws=NULL) {
    sum(cluster_loglik(theta, cd, method, K, draws))
}
