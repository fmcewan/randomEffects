#' Finding cluster modes and curvatures
#'
#' Finding the mode of each cluster's intercept and the
#' negative second derivative of its log-integrand (objective) at that mode
#'
#' @param theta Vector containing: 
#'   one regression coefficient per column, 
#'   followed by log(sigma), where sigma is the standard deviation of the random intercepts
#' @param cd List containing: 
#'   y: a numeric vector of binary responses,
#'   X: a finite numeric model matrix with one row per response,
#'   group: Holds integer identifiers (from 1 to J) matching patient to practice,
#'   and J: the number of clusters (int)
#' @param tol First derivative tolerance threshold, default = 1e-8
#' @param maxit Positive integer maximum number of Newton updates, default = 100
#'
#' @details
#' All cluster modes are updated simultaneously using Newton's method with analytic derivatives
#' Backtracking is applied separately to each cluster to ensure sufficient increase in objective
#' The input coefficients and random-intercept standard deviation are fixed
#' Convergence requires every absolute first derivative to be at most tol
#'
#' @return A list containing:
#' - u_hat: The mode for each practice, in cluster order
#' - H: The negative second derivative at each mode
#' - iterations: The number of Newton updates
#' - converged: TRUE if convergence is reached, throws error otherwise
#'
#' @examples
#' cd <- list(
#'   y = c(0, 1, 1, 1),
#'   X = cbind(1, c(-1, 0, 0.5, 1)),
#'   group = c(1, 1, 2, 2),
#'   J = 2
#' )
#' theta <- c(-0.5, 1, log(0.8))
#' fit <- cluster_modes(theta, cd)
#' fit$u_hat
#' fit$H
#'
#' @export
cluster_modes <- function(theta, cd, tol = 1e-8, maxit = 100) {
  
  # Checking that the expected inputs from cluster_data() exist
  if (!is.list(cd) || !all(c("y", "X", "group", "J") %in% names(cd))) {
    stop("cd must contain y, X, group and J.")
  }
  
  # Storing in new variables for simplicity
  y <- cd$y
  X <- cd$X
  group <- cd$group
  J <- cd$J
  
  # Checking more properties of input before running (see error messages)
  if (!is.matrix(X) || !is.numeric(X) || nrow(X) < 1 || ncol(X) < 1 || any(!is.finite(X))) {
    stop("cd$X must be a non-empty finite numeric matrix.")}
  if (!is.numeric(y) || length(y) != nrow(X) || anyNA(y) || any(!y %in% c(0, 1))) {
    stop("cd$y must contain one binary response per row of cd$X.")}
  if (!is.numeric(J) || length(J) != 1 || !is.finite(J) || J < 1 || J != floor(J)) {
    stop("cd$J must be a positive integer.")}
  if (!is.numeric(group) || length(group) != length(y) ||
      any(!is.finite(group)) || any(group != floor(group)) ||
      any(group < 1 | group > J) || length(unique(group)) != J) {
    stop("cd$group must use every integer cluster index from 1 to J.")}
  if (!is.numeric(theta) || length(theta) != ncol(X) + 1 ||
      any(!is.finite(theta))) {
    stop("theta must contain ncol(cd$X) coefficients followed by log(sigma).")}
  if (!is.numeric(tol) || length(tol) != 1 || !is.finite(tol) || tol <= 0 ||
      !is.numeric(maxit) || length(maxit) != 1 ||
      !is.finite(maxit) || maxit < 1 || maxit != floor(maxit)) {
    stop("tol must be positive and maxit must be a positive integer.")}
  
  # Extracting the regression coefficients and 1/sig^2
  beta <- theta[seq_len(ncol(X))]
  inv_sigma2 <- exp(-2 * theta[length(theta)])
  
  # Stop if inverse sigma2 computed is unusable 
  if (!is.finite(inv_sigma2) || inv_sigma2 <= 0) {
    stop("log(sigma) is outside the supported numerical range.")
  }
  
  # Vector of each patient's linear predictor eta (without the practice effect) 
  eta_fixed <- drop(X %*% beta)
  
  if (any(!is.finite(eta_fixed))) {
    stop("X %*% beta must be finite. Check parameter and covariate scales.")
  }
  
  # Defining some helper functions:
  
  # To sum patient contributions (per cluster)
  cluster_sum <- function(z) { drop(rowsum(z, group, reorder = TRUE)) }
  
  # Evaluating log(1 + exp(z)) without overflow (from chapter 10 notes)
  softplus <- function(z) { pmax(z, 0) + log1p(exp(-abs(z))) }
  
  # Evaluating each practice's log-integrand at proposed effects u (adapted idea from practical 7 h_j())
  log_integrand <- function(u) {
    # First, adding each practice effect to the predictors of its patients
    eta <- eta_fixed + u[group]
    
    log_prob <- ifelse(
      y == 1,
      -softplus(-eta), # If attending
      -softplus(eta)   # If not attending
    )
    # Penalizing too large (unlikely) effects, based on normal dist. assumption
    cluster_sum(log_prob) + dnorm(u, mean = 0, sd = exp(theta[length(theta)]), log = TRUE)
  }
  
  # For derivatives needed for Newthon's method
  derivatives <- function(u) {
    eta <- eta_fixed + u[group]
    
    # Computing prob of attendance and complement 
    p <- plogis(eta); q <- plogis(-eta) # Could have done q=1-p but risked rounding errors
    
    # Closed form derivatives as stated in assignment
    list(
      gradient = cluster_sum(ifelse(y == 1, q, -p)) - u * inv_sigma2,
      H = cluster_sum(p * q) + inv_sigma2
    )
  }
  
  u <- numeric(J)  # initialising practice effects as a vector of 0s of length J
  iterations <- 0
  
  repeat {
    d <- derivatives(u)
    
    # Checking derivatives/curvature for Newton update
    if (any(!is.finite(d$gradient)) ||
        any(!is.finite(d$H)) || any(d$H <= 0)) {
      stop("Non-finite derivatives or non-positive curvature encountered.")
    }
    
    # Can stop when derivatives are flat enough
    if (all(abs(d$gradient) <= tol)) {
      break
    }
    
    # Stop if max iterations reached
    if (iterations >= maxit) {
      stop("cluster_modes() did not converge; increase maxit or check inputs.")
    }
    
    active <- abs(d$gradient) > tol                  # Mask for non-optimised gradients
    step <- numeric(J)                               # Either 0 step update
    step[active] <- d$gradient[active] / d$H[active] # Or where needed store Newton step
    
    # Backtracking separately for each practice
    old_value <- log_integrand(u)
    accepted <- !active # Store non-active practices as already accepted 
    
    # Looping over backtracking attempts 
    for (bt in seq_len(60)) {
      candidate <- u + step # Trying step on all practices 
      candidate_value <- log_integrand(candidate) # Evaluating each practice for given effect
      sufficient_increase <- old_value + 1e-4 * d$gradient * step # Establish sufficient increase condition
      
      # Consider the step accepted if already accepted or met sufficient increase condition
      accepted <- accepted | (is.finite(candidate_value) & candidate_value >= sufficient_increase)
      
      # Stop when all steps accepted
      if (all(accepted)) {
        break
      }
      
      # Halve not accepted steps and keep going
      step[!accepted] <- step[!accepted] / 2
    }
    
    # Stop if loop failed 
    if (!all(accepted)) {
      stop("Newton backtracking failed to improve a cluster log-integrand.")
    }
    
    # Update effects
    u <- u + step
    iterations <- iterations + 1
  }
  
  # Return estimated practice effects, curvature, iterations and convergence flag
  list(u_hat = u, H = d$H, iterations = iterations, converged = TRUE)
}