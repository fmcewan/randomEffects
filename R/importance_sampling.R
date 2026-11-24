#' Standardised t draws for importance sampling
#' 
#' Returns a J x R matrix of draws from a standard t distribution, one row per cluster. Later 
#' shifted and scaled to give importance sampling draws of each clusters random effect. The same
#' draws are used everytime the likelihood is evaluated
#' 
#' @param cd List of clustered data. Note only cd$J, the number of clusters is used
#' @param R Number of draws per cluster. Must be at least 2 as importance sampling estimate is
#' an average over draws, so need 2 or more to measure MC variability.
#' @param df Degrees of freedom of t distribution. Small df means heavy tails, so the proposal
#' covers the tails of the target safely & the importance weightings stay stable.
#' 
#' @return A `cd$J` by `R` numeric matrix of draws, with attribute `"df"`.
#' @export
#' 
#' @examples
#' cd <- list(J = 5)
#' set.seed(1)
#' Z <- is_draws(cd, R = 100, df = 5)
#' dim(Z)
#' attr(Z, "df")
is_draws <- function(cd, R, df = 5) {
  
  # Check the inputs
  if (!is.list(cd) || is.null(cd$J)) stop("`cd` must be a list containing `J` ")
  if (!is.numeric(R) || length(R) != 1 || !is.finite(R) || R < 2 || R != round(R))
    stop("`R` must be a single whole number of at least 2.")
  if (!is.numeric(df) || length(df) != 1 || is.na(df) || df <= 0)
    stop("`df` must be a single positive number.")
  
  # One row per cluster, R columns of standard t draws
  Z <- matrix(stats::rt(cd$J * R, df = df), nrow = cd$J, ncol = R)
  
  # Record df on the matrix, then return it
  attr(Z, "df") <- df
  Z
}

# Helper functions for importance_sampling
# log of the curve height, h_i(b), at a J x R matrix of spots B (row i = practice i)

# This function was AI assisted/adapted: helped suggest most efficient way to create matrix
log_integrand_matrix <- function(B, beta, sigma, cd) {
  # stable log(1 + e^z), as suggested by Practical 5 
  softplus <- function(z) pmax(z, 0) + log1p(exp(-abs(z)))  
  # n x R: each patient gets its practice's spots
  eta <- drop(cd$X %*% beta) + B[cd$group,,drop = FALSE]    
  # J x R: sum over each practice's patients
  loglik <- unname(rowsum(cd$y * eta - softplus(eta), cd$group, reorder = TRUE))  
  loglik + stats::dnorm(B, 0, sigma, log = TRUE)
}

# This function was AI assisted: improved my version of how assignment instructed to implement log-sum-exp
# log(rowSums(exp(A))) without underflow: factor out each row's largest value
row_logsumexp <- function(A) {
  amax <- A[cbind(seq_len(nrow(A)), max.col(A, ties.method = "first"))]
  amax + log(rowSums(exp(A - amax)))
}

#' Importance Sampling estimate of each clusters log-likelihood
#'
#' importance_sampling() was adapted from code from Practical 7
#' 
#' Estimates log L_i = log integral with respect to b (f(y_i|b)*phi(b;0, sigma))
#' for each cluster. The proposal for cluster i is a t distribution centred at
#' the mode of that cluster's integrand and scaled by its curvature,
#' b = uhat_i + T/root(H_i). L_i is then estimated by the mean
#' importance weight f(b)/q(b), with everything computed on the log scale.
#'
#' @param theta Numeric vector (beta, logsigma).
#' @param cd List with `y` (0/1 response), `X` (n x p design matrix), `group`
#'   (cluster index 1..J) and `J` (number of clusters).
#' @param draws J x R matrix from [is_draws()]. Reuse the same matrix for every
#'   `theta`, so that the estimate stays smooth in `theta` and can be optimised.
#'
#' @return Numeric vector of length J of estimated logL_i, with
#'   attribute `"mc_se"` giving the Monte Carlo standard error of each.
#'
#' @seealso [is_draws()], [cluster_modes()]
#' @keywords internal
#'
#' @seealso [is_draws()], [cluster_modes()]
#' @noRd
importance_sampling <- function(theta, cd, draws) {
  #Check the draws are usable
  if (is.null(draws) || !is.matrix(draws) || nrow(draws) != cd$J || is.null(attr(draws, "df")))
    stop("method 'is' needs `draws` from is_draws(cd, R, df), with one row per cluster.")
  
  # Unpack theta = (beta, log sigma)
  p <- ncol(cd$X)
  beta <- theta[seq_len(p)]
  sigma <- exp(theta[p + 1])
  R <- ncol(draws)
  
  # centre and width of each practices spike
  md <- cluster_modes(theta, cd)
  centre <- unname(md$u_hat)
  width <- 1 / sqrt(unname(md$H))
  
  #Place the spots (J x R)
  B <- centre + width * draws
  
  # Log height of the curve at each spot
  log_f <- log_integrand_matrix(B, beta, sigma, cd)
  
  # Log proposal density at each spot
  log_q <- stats::dt(draws, df = attr(draws, "df"), log = TRUE) - log(width)
  
  # Log weights
  log_w <- log_f - log_q
  
  # Log of the average weight, per practice
  log_L <- row_logsumexp(log_w) - log(R)
  
  # Monte Carlo se of each log L_i (from how much the weights vary)
  w_rel <- exp(log_w - log_L) # weights relative to their average (these average to 1)
  attr(log_L, "mc_se") <- sqrt(rowSums((w_rel - 1)^2) / (R - 1)) / sqrt(R)
  log_L
}