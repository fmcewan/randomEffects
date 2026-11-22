# negative_loglik
negative_loglik <- function(theta, cd, lik_method,  K, draws) {
  -marginal_loglik(theta, cd, lik_method, K, draws)
}

#' Fit a random intercept logistic regression
#'
#' This method maximises the marginal log-likelihood over theta=(beta,log(sigma))
#' using BFGS. Standard errors are generated from a numerical Hessian calculated
#' using numDeriv::hessian.
#'
#' @param cd Data prepared by cluster_data.
#' @param method Chosen integration method ("gh", "laplace", "agh" or "is").
#' @param K Number of quadrature nodes (used by "gh" and "agh").
#' @param draws Importance-sampling draws from is_draws.
#' @param start Starting values for theta. 
#'              By default, a pooled logistic regression for beta 
#'              and log sigma = 0.
#'
#' @return A list with par (estimates of theta), se (standard
#'   errors, which are NA if the Hessian cannot be inverted), loglik
#'   (the maximised log-likelihood) and convergence (code from
#'   optim, 0 means success).
#'
#' @importFrom stats optim glm coef binomial
#' @export 
#' @examples
#' data <- rint_sim(J = 50, m = 5, beta = c(-0.5, 1), sigma = log(1))
#' cd <- cluster_data(data$y, data$X, factor(data$group))
#' fit_rint(cd, method = "laplace")
fit_rint <- function(cd, method, K=NULL, draws=NULL, start=NULL) {
  
  # If no start is supplied, guess through logistic regression
  if (is.null(start)) {
    beta_start <- unname(coef(glm(cd$y ~ cd$X - 1, family = binomial)))
    start <- c(beta_start, 0)
  }
  
  #AI suggested fix
  # objective for optim(): treat points where the likelihood cannot be evaluated as infinitely bad,
  #so BFGS backtracks instead of stopping
  safe_nll <- function(theta, ...) {
    value <- tryCatch(negative_loglik(theta, ...), error = function(e) Inf)
    if (is.finite(value)) value else Inf
  }
  
  # Define optimisers and Hessian
  optimiser <- optim(start, safe_nll, cd=cd, lik_method=method, K=K, draws=draws, method="BFGS")
  hessian <- numDeriv::hessian(negative_loglik, optimiser$par, cd=cd, lik_method=method, K=K, draws=draws)

  # Find standard errors
  standard_errors <- tryCatch(sqrt(diag(solve(hessian))), error = function(e) rep(NA, length(optimiser$par)))

  return (list(par=optimiser$par, se=standard_errors, loglik=-optimiser$value, convergence=optimiser$convergence))  

}

