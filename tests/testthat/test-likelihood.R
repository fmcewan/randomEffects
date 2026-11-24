# Fake cluster_logliks
fake_loglik <- function(...) c(-1, -2, -3)
fake_loglik_with_se <- function(...) {
  structure(c(-1, -2, -3), mc_se = rep(0.1, 3))
}

# Unit tests 
test_that("Invalid method types are rejected in cluster_loglik", {
  expect_error(cluster_loglik(theta, cd, method=NA, K=1, draws=NULL), "Invalid method name.")
  expect_error(cluster_loglik(theta, cd, method=123, K=1, draws=NULL), "Invalid method name.")
  expect_error(cluster_loglik(theta, cd, method=c("gh", "agh"), K=1, draws=NULL), "Invalid method name.")
  expect_error(cluster_loglik(theta, cd, method=NULL, K=1, draws=NULL), "Invalid method name.")
})

test_that("Unknown method is rejected in cluster_loglik", {
  expect_error(cluster_loglik(theta, cd, method="unknown", K=1, draws=NULL), "Method not yet implemented.")
})

test_that("The cluster log-likelihoods are summed by marginal_loglik", {
  testthat::local_mocked_bindings(cluster_loglik = fake_loglik)  

  expect_equal(marginal_loglik(theta, cd, method="laplace"), -6)
})

test_that("marginal_loglik returns a single number", {
  testthat::local_mocked_bindings(cluster_loglik = fake_loglik_with_se)

  out <- marginal_loglik(theta, cd, method="is")
  
  expect_length(out, 1)
  expect_null(attributes(out))
})

# Integration tests
# TODO: Integration tests for gh

test_that("cluster_loglik dispatches to the real laplace method", {
  out <- cluster_loglik(theta, cd, method = "laplace")

  expect_length(out, cd$J)
  expect_true(all(is.finite(out)))
})

test_that("cluster_loglik dispatches to the real agh method", {
  out <- cluster_loglik(theta, cd, method = "agh", K = 7)

  expect_length(out, cd$J)
  expect_true(all(is.finite(out)))
})

test_that("marginal_loglik equals the sum of the real cluster log-likelihoods", {
  expect_equal(marginal_loglik(theta, cd, method = "laplace"), 
    sum(cluster_loglik(theta, cd, method = "laplace")))
  expect_equal(marginal_loglik(theta, cd, method = "agh", K = 7), 
    sum(cluster_loglik(theta, cd, method = "agh", K = 7)))
})

test_that("agh with lots of nodes agrees with integrate()", {
  agh <- cluster_loglik(theta, cd, method = "agh", K = 30)
  expect_equal(unname(agh), true_loglik(theta, cd), tolerance = 1e-6)
})

test_that("laplace agrees with integrate()", {
  laplace <- cluster_loglik(theta, cd, method = "laplace")
  expect_equal(unname(laplace), true_loglik(theta, cd), tolerance = 0.05)
})

test_that("agh with one node equals laplace", {
  expect_equal(cluster_loglik(theta, cd, method = "agh", K = 1), cluster_loglik(theta, cd, method = "laplace"))
})

test_that("methods that need K or draws give an error without them", {
  expect_error(cluster_loglik(theta, cd, method = "agh", K = NULL), "positive integer")
  expect_error(cluster_loglik(theta, cd, method = "is", draws = NULL), "needs")
})

test_that("importance sampling is within a few MC standard errors of integrate()", {
  set.seed(101)
  draws <- is_draws(cd, R = 1000)

  estimate <- cluster_loglik(theta, cd, method = "is", draws = draws)
  truth <- true_loglik(theta, cd)

  expect_false(is.null(attr(estimate, "mc_se")))
  expect_true(all(abs(as.numeric(estimate) - truth) < 4 * attr(estimate, "mc_se")))
})
