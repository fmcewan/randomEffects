# Tests for cluster_modes()

# Test data
cm_cd <- list(
    y = c(0, 0, 0, 0,  # Practice 1 has all-zero responses
          1, 1, 1, 1,  # Practice 2 has all-one responses
          0, 1, 0, 1), # And practice 3 has a mixture
    X = cbind(1, rep(c(-1, -0.5, 0.5, 1), 3)), group = rep(1:3, each = 4), J = 3)
# Parameter theta to run cluster_modes
cm_theta <- c(-0.5, 1, log(0.8))
# Extracting these from theta to avoid repetition
beta <- cm_theta[1:2]
sigma <- exp(cm_theta[3])

# TEST 1) Checking if gradient met threshold, and convergence
test_that("cluster modes have zero gradient", {
    fit <- cluster_modes(cm_theta, cm_cd)
    eta <- drop(cm_cd$X %*% beta) + fit$u_hat[cm_cd$group]
    p <- plogis(eta)
    gradient <- drop(rowsum(cm_cd$y - p, cm_cd$group, reorder = TRUE)) - fit$u_hat / sigma^2

    # First checking converged, and appropriateness of outputs u and H
    expect_true(fit$converged)
    expect_length(fit$u_hat, cm_cd$J)
    expect_length(fit$H, cm_cd$J)
    expect_true(all(is.finite(fit$u_hat)))
    expect_true(all(is.finite(fit$H)))
    expect_true(all(fit$H > 0))

    # Checking derivative is small enough
    expect_lt(max(abs(gradient)), 1e-7)
})


# TEST 2) Verifying with independent calculation
# Maximising with optimize() then computing Hessian (adapted from Practical 7)
test_that("modes and curvatures agree with independent calculations", {
    fit <- cluster_modes(cm_theta, cm_cd)

    # Looping over practices (just for testing)
    for (j in seq_len(cm_cd$J)) {
        patients <- cm_cd$group == j

        # Log-integrand for this practice only (adapted h_j() from practical 7)
        h_j <- function(u) {
            eta <- drop(cm_cd$X[patients, , drop = FALSE] %*% beta) + u
            sum(dbinom(cm_cd$y[patients], size = 1, prob = plogis(eta), log = TRUE)) + 
            dnorm(u, mean = 0, sd = sigma, log = TRUE)
        }
        # Using optimize and numDeriv to maximise h_j and compute its curvature
        reference <- optimize(h_j, interval = c(-10, 10), maximum = TRUE, tol = 1e-10)
        reference_H <- -numDeriv::hessian(h_j, reference$maximum)[1, 1]

        # Checking if u,H found by cluster_modes are close enough to the one found using optimize
        expect_equal(unname(fit$u_hat[j]), reference$maximum, tolerance = 1e-5)
        expect_equal(unname(fit$H[j]), reference_H, tolerance = 1e-5)
    }
})


# TEST 3) Check "extreme-valued" practices with larger variation.
test_that("extreme response patterns work with larger sigma", {
    # Defining new theta (larger variation)
    theta <- c(-0.5, 1, log(5))
    fit <- cluster_modes(theta, cm_cd)

    # Standard checks as in Test 1
    expect_true(fit$converged)
    expect_true(all(is.finite(fit$u_hat)))
    expect_true(all(is.finite(fit$H)))
    expect_true(all(fit$H > 0))

    expect_lt(fit$u_hat[1], 0) # All-zero practice should have a negative effect
    expect_gt(fit$u_hat[2], 0) # All-one practice should have a positive effect
})


# TEST 4) Testing errors from invalid inputs 
test_that("invalid inputs produce informative errors", {

  # Required fields missing (empty cd)
  expect_error(cluster_modes(cm_theta, list()), regexp = "cd must contain")

  # Non binary response
  bad_cd <- cm_cd
  bad_cd$y[1] <- 2 # replacing with one non 0-1 value
  expect_error(cluster_modes(cm_theta, bad_cd), regexp = "binary response")

  # Inappropriate theta provided (no log(sigma))
  expect_error(
    cluster_modes(cm_theta[-3], cm_cd), # removing log(sigma) from theta before running 
    regexp = "theta must contain"
  )

  # Negative tolerance
  expect_error(cluster_modes(cm_theta, cm_cd, tol = -1), regexp = "tol must be positive")
})