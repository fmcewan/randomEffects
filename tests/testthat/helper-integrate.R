# Reference loglik function 
integrate_cluster <- function(u, y, eta0, sigma) {
  sapply(u, function(ui) {
    eta <- eta0 + ui
    exp(sum(y * eta - log1p(exp(eta))) + dnorm(ui, mean = 0, sd = sigma, log = TRUE))
  })
}

true_loglik <- function(theta, cd) {
  beta  <- theta[seq_len(ncol(cd$X))]
  sigma <- exp(theta[length(theta)])
  eta0  <- drop(cd$X %*% beta)

  result <- numeric(cd$J)

  for (j in seq_len(cd$J)) {
    in_cluster <- cd$group == j

    integral <- integrate(integrate_cluster, -Inf, Inf,
                          y = cd$y[in_cluster],
                          eta0 = eta0[in_cluster],
                          sigma = sigma)

    result[j] <- log(integral$value)
  }

  result
}
