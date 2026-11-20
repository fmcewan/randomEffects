devtools::load_all()     # gives access to internal functions like importance_sampling()

# One data set (same setup as B1) 
set.seed(4113)  
beta_true <- c(-0.5, 1)
sigma_true <- 1.5
theta0 <- c(beta_true, log(sigma_true))

sim  <- rint_sim(50, 20, beta_true, log(sigma_true))   
cd20 <- cluster_data(sim$y, sim$X, factor(sim$group))  # cluster_data needs a factor

# Truth
truth <- sum(agh_loglik(theta0, cd20, K = 50))

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

