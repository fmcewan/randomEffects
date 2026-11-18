# Part B1, B2, B3

# B1



# B2



# B3: Functions and settings for computational cost comparisons

# Default settings (to be able to pass this as one variable later)
b3_settings <- list(
  cluster_counts = c(50, 200, 800),
  m = 5,                          # Patients per practice
  beta = c(-0.5, 1), sigma = 1.5, # Theta
  K_gh = 10, K_agh = 7,           # Nodes
  R_is = 100, df_is = 4,          # For importance sampling
  min_seconds = 0.2,              # For stable timings
  rounds = 3                      # Num rounds to be timed
)

# Timed evaluations with system.time()
b3_time_evaluation <- function(f, settings = b3_settings) {
  min_seconds <- settings$min_seconds
  rounds <- settings$rounds

  invisible(f())    # Checking the call once without returning or timing
  repetitions <- 1  # Starting with 1 call
  # Defining helper to measure elapsed time
  measure_elapsed <- function(n) {
    system.time({
      for (i in seq_len(n)) {
        invisible(f())
      }
    })[["elapsed"]]
  }

  # Repeating until proper batch size is found, then timing whole batch
  repeat {
    elapsed <- measure_elapsed(repetitions)
    if (elapsed >= min_seconds) {
      break # break if long enough
    }
    repetitions <- repetitions * 2 # Otherwise increase num of calls and retry
  }

  # And repeating batches to show the variability in the timings
  seconds_per_call <- numeric(rounds)
  for (r in seq_len(rounds)) {
    elapsed <- measure_elapsed(repetitions)
    seconds_per_call[r] <- elapsed / repetitions # Finding time per repetition
  }

  # Individual batch timings and summaries
  list(
    median_seconds = median(seconds_per_call),
    min_seconds = min(seconds_per_call),
    max_seconds = max(seconds_per_call),
    calls_per_batch = repetitions,
    seconds_per_call = seconds_per_call
  )
}

# Calling optimize() separately for every practice (just for timing purposes)
b3_separate_modes <- function(theta, cd) {

  beta <- theta[seq_len(ncol(cd$X))];  sigma <- exp(theta[length(theta)])
  eta_fixed <- drop(cd$X %*% beta)
  u_hat <- H <- numeric(cd$J) # Both 0 vectors of len J

  softplus <- function(z) {
    pmax(z, 0) + log1p(exp(-abs(z)))
  }

  for (j in seq_len(cd$J)) {

    patients <- cd$group == j
    y_j <- cd$y[patients]
    eta_j <- eta_fixed[patients]

    h_j <- function(u) {
      eta <- eta_j + u
      log_prob <- ifelse(y_j == 1, -softplus(-eta), -softplus(eta))
      sum(log_prob) + dnorm(u, mean = 0, sd = sigma, log = TRUE)
    }

    # Bound to find modes based on the derivative equation (+ margin of 1)
    bound <- length(y_j) * sigma^2 + 1
    opt <- optimize(h_j, interval = c(-bound, bound),
                    maximum = TRUE, tol = 1e-8)
    u_hat[j] <- opt$maximum

    # Using analytic curvature to reduce calculations
    eta <- eta_j + u_hat[j]
    # Storing the curvature for practice j's mode
    H[j] <- sum(plogis(eta) * plogis(-eta)) + 1 / sigma^2
  }

  list(
    u_hat = setNames(u_hat, seq_len(cd$J)),
    H = setNames(H, seq_len(cd$J))
  )
}

# Both comparisons needed for B3
# (+ returns results to display directly in the report)
run_b3_cost <- function(settings = b3_settings) {

  theta <- c(settings$beta, log(settings$sigma))

  # Preparing datasets and fixed draws before starting timing
  datasets <- lapply(settings$cluster_counts, function(J) {
    rint_sim(J = J, m = settings$m, beta = settings$beta,
             sigma = settings$sigma)
  })
  draws <- lapply(datasets, function(cd) {
    is_draws(cd, R = settings$R_is, df = settings$df_is)
  })

  likelihood_rows <- mode_rows <- list()

  for (i in seq_along(datasets)) {

    J <- settings$cluster_counts[i]; cd <- datasets[[i]]

    # Comparing evaluations for all methods
    for (method in c("gh", "laplace", "agh", "is")) {

      # K only for the quadrature methods
      K <- if (method == "gh") settings$K_gh else settings$K_agh
      method_draws <- if (method == "is") draws[[i]] else NULL

      timing <- b3_time_evaluation(
        function() {
                    marginal_loglik(theta, cd, method = method,
                                    K = K, draws = method_draws)},
        settings = settings
      )

      likelihood_rows[[length(likelihood_rows) + 1]] <- data.frame(
        J = J,
        method = method,
        median_seconds = timing$median_seconds,
        min_seconds = timing$min_seconds,
        max_seconds = timing$max_seconds,
        calls_per_batch = timing$calls_per_batch
      )
    }

    # Confirm that both mode finders produce comparable answers.
    vectorised <- cluster_modes(theta, cd)
    separate <- b3_separate_modes(theta, cd)

    mode_difference <- max(abs(vectorised$u_hat - separate$u_hat))
    curvature_difference <- max(abs(vectorised$H - separate$H))

    if (!is.finite(mode_difference) || !is.finite(curvature_difference) ||
          mode_difference >= 1e-5 || curvature_difference >= 1e-5) {
      stop("Mode-finding implementations disagree for J = ", J)
    }

    vectorised_time <- b3_time_evaluation(function() cluster_modes(theta, cd),
                                          settings = settings)
    separate_time <- b3_time_evaluation(function() b3_separate_modes(theta, cd),
                                        settings = settings)

    mode_rows[[i]] <- data.frame(
      J = J,
      vectorised_seconds = vectorised_time$median_seconds,
      separate_seconds = separate_time$median_seconds,
      # speed-up > 1 means the vectorised version is faster
      speed_up = separate_time$median_seconds / vectorised_time$median_seconds,
      max_mode_difference = mode_difference,
      max_curvature_difference = curvature_difference
    )
  }

  # Return tables for the report
  list(
    settings = settings,
    likelihood_timings = do.call(rbind, likelihood_rows),
    mode_timings = do.call(rbind, mode_rows)
  )
}