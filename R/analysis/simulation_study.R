# Part C1, C2, C3, C4

specifications <- data.frame(
  method = c("laplace", "gh", "agh", "is"),
  K = c(NA, 10, 7, NA),
  R = c(NA, NA, NA, 100)
)

fit_method <- function(cd, method, truth, K = NULL, draws = NULL) {

  # Time how long it takes fit_rint to run
  start <- Sys.time()
  fit <- tryCatch(
    fit_rint(cd, method = method, K = K, draws = draws),
    error = function(e) NULL
  )
  time <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  if (is.null(fit)) {
    theta <- rep(NA, 3)
    standard_error_theta <- rep(NA, 3)
    converged <- FALSE
  } else {
    theta <- fit$par
    standard_error_theta <- fit$se
    converged <- fit$convergence == 0 && all(is.finite(standard_error_theta))
  }

  b1_estimate <- theta[2]
  b1_standard_error <- standard_error_theta[2]
  b1_lower <- b1_estimate - 1.96 * b1_standard_error
  b1_upper <- b1_estimate + 1.96 * b1_standard_error

  sigma_estimate <- exp(theta[3])
  sigma_standard_error <- sigma_estimate * standard_error_theta[3]
  sigma_lower <- exp(theta[3] - 1.96 * standard_error_theta[3])
  sigma_upper <- exp(theta[3] + 1.96 * standard_error_theta[3])

  data.frame(
    method = method,
    param = c("beta1", "sigma"),
    true = unname(truth),
    est = c(b1_estimate, sigma_estimate),
    se = c(b1_standard_error, sigma_standard_error),
    lower = c(b1_lower, sigma_lower),
    upper = c(b1_upper, sigma_upper),
    converged = converged,
    time = time
  )
}
