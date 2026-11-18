devtools::load_all()     # gives access to internal functions like importance_sampling()

# One data set (same setup as B1) 
set.seed(4113)  
beta_true <- c(-0.5, 1)
sigma_true <- 1.5
theta0 <- c(beta_true, log(sigma_true))
#! TEMPORARY: swap for rint_sim() + cluster_data() when merged !
cd20 <- sim_cd(50, 20, beta_true, sigma_true)  

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
    predicted = -0.5 * mean(reps[, "mc_se"]^2)     # Jensen: bias approx -SE^2 / 2
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
       y = "Bias of IS log-likelihood",
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



# Plotting currently fails: cluster_modes() does not converge for all sigma
set.seed(4114)
sigma_grid <- seq(1, 2, length.out = 60)
Z_fixed <- is_draws(cd20, R = 100)                          # one draw set, reused everywhere

profile <- data.frame(
  sigma = rep(sigma_grid, 2),
  ell = c(sapply(sigma_grid, function(s)                          # Same draws at every sigma
      sum(importance_sampling(c(beta_true, log(s)), cd20, Z_fixed))),
      sapply(sigma_grid, function(s)                          #New draws at every sigma
      sum(importance_sampling(c(beta_true, log(s)), cd20, is_draws(cd20, R = 100))))
  ),
  draws = rep(c("Same draws (common random numbers)", "Fresh draws each evaluation"),
              each = length(sigma_grid))
)

p_crn <- ggplot(profile, aes(x = sigma, y = ell, colour = draws)) +
  geom_line() +
  labs(x = expression(sigma ~ "(beta fixed at the truth)"),
       y = "IS log-likelihood (R = 100)", colour = NULL,
       title = "Fixed draws give a smooth objective. Fresh draws give a jagged one") +
  theme_bw() + theme(legend.position = "bottom")
p_crn

