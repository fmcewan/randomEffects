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

  # Initiaise necessary variables; if they don't exist, NA 
  if (is.null(fit)) {
    theta <- rep(NA, 3)
    standard_error_theta <- rep(NA, 3)
    converged <- FALSE
  } else {
    theta <- fit$par
    standard_error_theta <- fit$se
    converged <- fit$convergence == 0 && all(is.finite(standard_error_theta))
  }

  # Calculations 
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

simulate_replicate <- function(replicate, scenario, J, beta, seed) {
  
  # Set unique seed for each replicate 
  set.seed(seed+100000*scenario$id+replicate)
  
  # Define necessary variables 
  simulation <- rint_sim(J, scenario$m, beta, scenario$sigma)
  cd <- cluster_data(simulation$y, simulation$X, simulation$group)
  R = specifications$R[specifications$method == "is"] 
  draws <- is_draws(cd, R=R)
  truth <- c(beta[2], scenario$sigma)
  
  # Fit each method in the specification
  fits <- vector("list", nrow(specifications))
  for (i in seq_len(nrow(specifications))) {
    fits[[i]] <- fit_method(cd, specifications$method[i], truth, K=specifications$K[i], draws=draws)
  }
  
  # Stack the results into one table 
  results <- do.call(rbind, fits)
  
  # Put the scenaro deatils in front of every row 
  cbind(rep = replicate, m=scenario$m, sigma=scenario$sigma, results)
}

run_study <- function(number_of_replicates, J = 50, beta = c(-0.5, 1), seed = 2026) {
  
  # Create scenarios
  scenarios <- expand.grid(m = c(2, 5, 20), sigma = c(0.5, 1.5))
  scenarios$id <- seq_len(nrow(scenarios))
  
  # Create results vector
  results <- vector("list", nrow(scenarios) * number_of_replicates)
  k <- 0 
  
  # Simulate replicates for all results
  for (s in seq_len(nrow(scenarios))) {
    for (r in seq_len(number_of_replicates)) {
      k <- k + 1
      results[[k]] <- simulate_replicate(r, scenarios[s, ], J, beta, seed)
    }
  }

  do.call(rbind, results)

}
