
#B1--------
calculate_errors <- function(sim_data, theta, K) {
  
  # Use adaptive Gauss-Hermite with 50 nodes as the reference
  agh_50 <- marginal_loglik(theta, sim_data, "agh", K = 50)
  
  # Calculate both errors for every value of K
  integration_errors <- vapply(
    K,
    function(k) {
      c(
        gh_error = abs(marginal_loglik(theta, sim_data, "gh", k) - agh_50),
        agh_error = abs(marginal_loglik(theta, sim_data, "agh", k) - agh_50)
      )
    },
    numeric(2)
  )
  
  # Calculate Laplace error once, since it does not depend on K
  laplace_error <- abs(
    marginal_loglik(theta, sim_data, "laplace") - agh_50
  )
  
  # Combine everything into a data frame
  data.frame(
    K = K,
    gh_error = integration_errors["gh_error", ],
    agh_error = integration_errors["agh_error", ],
    laplace_error = rep(laplace_error, length(K))
  )
}

plot_errors <- function(errors, m) {
  
  ggplot(errors, aes(x = K)) +
    geom_path(aes(y = gh_error, colour = "Gauss-Hermite")) +
    geom_path(aes(y = agh_error, colour = "Adaptive Gauss-Hermite")) +
    geom_path(aes(y = laplace_error, colour = "Laplace")) +
    scale_colour_manual(
      values = c(
        "Gauss-Hermite" = "green",
        "Adaptive Gauss-Hermite" = "red",
        "Laplace" = "blue"
      )
    ) +
    scale_y_log10() +
    labs(
      x = "Number of integration nodes (K)",
      y = "Absolute error",
      colour = "Method",
      title = paste("Integration errors (m =", m, ")")
    ) 
}

# Simulation settings
set.seed(4113)
beta <- c(1, 0.2)
theta <- c(beta, log(1.5))
K <- 1:30

# Generate datasets
sim_data_20 <- rint_sim(50, 20, beta, log(1.5))
sim_data_2 <- rint_sim(50, 2, beta, log(1.5))

# Calculate errors
errors_20 <- calculate_errors(sim_data_20, theta, K)
errors_2 <- calculate_errors(sim_data_2, theta, K)

# Generate plots
p20 <- plot_errors(errors_20, 20)
p2 <- plot_errors(errors_2, 2)

# Display both plots
print(p20)
print(p2)




# B2----------
devtools::load_all()     # gives access to internal functions like importance_sampling()

# One data set (same setup as B1) 
set.seed(4113)  
beta_true <- c(-0.5, 1)
sigma_true <- 1.5
theta0 <- c(beta_true, log(sigma_true))

sim  <- rint_sim(50, 20, beta_true, log(sigma_true))   
cd20 <- cluster_data(sim$y, sim$X, factor(sim$group))  # cluster_data needs a factor

# Truth
truth <- marginal_loglik(theta0, cd20, method = "agh", K = 50)

#Repeat IS 200 times for each R
Rs <- c(10, 30, 100, 300, 1000)
n_rep <- 200

one_run <- function(R) {
  est <- importance_sampling(theta0, cd20, is_draws(cd20, R = R))
  c(ell = sum(est),                              # total IS loglikelihood
    mc_se = sqrt(sum(attr(est, "mc_se")^2)))        # package SE of that total
}

results <- do.call(rbind, lapply(Rs, function(R) {
  reps <- t(replicate(n_rep, one_run(R)))           # 200 rows: ell, mc_se
  data.frame(
    R = R,
    bias = mean(reps[, "ell"]) - truth,
    bias_mcse = sd(reps[, "ell"]) / sqrt(n_rep),   # uncertainty in the bias itself
    sd = sd(reps[, "ell"]),
    mean_mc_se = mean(reps[, "mc_se"]),
    predicted = -0.5 * mean(reps[, "mc_se"]^2)     # Jensen: bias approx -SE^2 / 2 (noted from lectures)
  )
}))
results


library(ggplot2)

p_bias <- ggplot(results, aes(x = R, y = bias)) +
  geom_hline(yintercept = 0, colour = "grey60") +           # zero bias reference
  geom_pointrange(aes(ymin = bias - 2 * bias_mcse,           # +-2 MCSE error bars
                      ymax = bias + 2 * bias_mcse)) +
  geom_line(aes(y = predicted), linetype = "dashed", colour = "steelblue") +  # Jensen prediction
  scale_x_log10(breaks = Rs) +
  labs(x = "Draws per cluster, R (log scale)",
       y = "Bias of IS log likelihood",
       title = "Bias of the IS log likelihood estimate (200 repetitions per R)",
       caption = "Error bars: +-2 Monte Carlo SEs. Dashed: Jensen prediction -SE^2/2.") +
  theme_bw()
p_bias

spread <- data.frame(
  R = rep(results$R, 2),
  value = c(results$sd, results$mean_mc_se),
  measure = rep(c("SD across 200 repetitions", "Mean reported MC SE"), each = nrow(results))
)

p_sd <- ggplot(spread, aes(x = R, y = value, colour = measure)) +
  geom_point() + geom_line() +
  scale_x_log10(breaks = Rs) + scale_y_log10() +
  labs(x = "Draws per cluster, R (log scale)", y = "Standard deviation (log scale)",
       colour = NULL,
       title = "Spread of the IS log-likelihood falls like 1/sqrt(R)") +
  theme_bw() + theme(legend.position = "bottom")
p_sd



sigma_grid <- seq(1, 2, length.out = 60)
Z_fixed <- is_draws(cd20, R = 100)                          # one draw set, reused everywhere

profile_df <- data.frame(
  sigma = rep(sigma_grid, 2),
  ell = c(sapply(sigma_grid, function(s)                          # Same draws at every sigma
    sum(importance_sampling(c(beta_true, log(s)), cd20, Z_fixed))),
    sapply(sigma_grid, function(s)                          #New draws at every sigma
      sum(importance_sampling(c(beta_true, log(s)), cd20, is_draws(cd20, R = 100))))
  ),
  draws = rep(c("Same draws (common random numbers)", "Fresh draws each evaluation"),
              each = length(sigma_grid))
)

p_crn <- ggplot(profile_df, aes(x = sigma, y = ell, colour = draws)) +
  geom_line() +
  labs(x = expression(sigma ~ "(beta fixed at the truth)"),
       y = "IS log likelihood (R = 100)", colour = NULL,
       title = "Fixed draws: smooth objective; fresh draws: jagged") +
  theme_bw() + theme(legend.position = "bottom")
p_crn



# Starting values as in fit_rint(): pooled logistic regression, log sigma = 0
start <- c(unname(coef(glm(cd20$y ~ cd20$X - 1, family = binomial))), 0)

# Incorrect- new draws every time the optimiser evaluates the likelihood
nll_fresh <- function(theta) -sum(importance_sampling(theta, cd20, is_draws(cd20, R = 100)))

#function was AI assisted - improved original to return data frame per fit and record conditionMessage(e)
# so that the reason each fresh draws fit failed is kept
fresh_fits <- do.call(rbind, lapply(1:5, function(k) {
  tryCatch({o <- optim(start, nll_fresh, method = "BFGS")
  se <- tryCatch(sqrt(diag(solve(numDeriv::hessian(nll_fresh, o$par)))),
                 error = function(e) rep(NA, 3))
  data.frame(fit = k, beta0 = o$par[1], beta1 = o$par[2], sigma = exp(o$par[3]),
             convergence = o$convergence, evals = unname(o$counts["function"]),
             se_beta1 = se[2], se_logsigma = se[3], error = NA)
  }, error = function(e) {
    data.frame(fit = k, beta0 = NA, beta1 = NA, sigma = NA,
               convergence = NA, evals = NA,
               se_beta1 = NA, se_logsigma = NA,
               error = conditionMessage(e))       # why it failed
  })
}))

fixed_fit <- fit_rint(cd20, method = "is", draws = is_draws(cd20, R = 100))

fixed_row <- data.frame(fit = "fixed", beta0 = fixed_fit$par[1], beta1 = fixed_fit$par[2],
                        sigma = exp(fixed_fit$par[3]), convergence = fixed_fit$convergence,
                        evals = NA, se_beta1 = fixed_fit$se[2],
                        se_logsigma = fixed_fit$se[3], error = NA)

rbind(fresh_fits, fixed_row)



#B3-------------
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
             sigma = log(settings$sigma))
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




