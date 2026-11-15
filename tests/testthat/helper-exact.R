#! THIS IS A STAND-IN UNTIL 'rint_sim()' EXISTS !

# Simulate a small clustered data set
sim_cd <- function(J, m, beta, sigma) {
  group <- rep(seq_len(J), each = m)
  x <- stats::rnorm(J * m)
  b <- stats::rnorm(J, 0, sigma)
  y <- stats::rbinom(J * m, 1, stats::plogis(beta[1] + beta[2] * x + b[group]))
  list(y = y, X = cbind(1, x), group = group, J = J)
}

# ~Exact logL_i for every cluster by integrate(), as a reference
exact_loglik <- function(theta, cd) {
  p <- ncol(cd$X); beta <- theta[1:p]; sigma <- exp(theta[p + 1])
  md <- cluster_modes(theta, cd)
  sapply(seq_len(cd$J), function(i) {
    rows <- cd$group == i
    lin <- drop(cd$X[rows, , drop = FALSE] %*% beta)
    h_i <- function(b) {
      eta <- outer(lin, b, "+")
      colSums(cd$y[rows] * eta - log1p(exp(eta))) + stats::dnorm(b, 0, sigma, log = TRUE)
    }
    mode <- unname(md$u_hat[i]); s <- 1 / sqrt(unname(md$H[i]))
    h0 <- h_i(mode)
    I <- stats::integrate(function(b) exp(h_i(b) - h0), mode - 30 * s,
                           mode + 30 * s, rel.tol = 1e-12)$value
    h0 + log(I)
  })
}