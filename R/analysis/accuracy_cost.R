# Part B1, B2, B3

#B1
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