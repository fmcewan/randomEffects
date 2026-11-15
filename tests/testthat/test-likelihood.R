# Fake cluster_logliks
fake_loglik <- function(...) c(-1, -2, -3)
fake_loglik_with_se <- function(...) {
  structure(c(-1, -2, -3), mc_se = rep(0.1, 3))
}

# TODO: Add integration tests once the real methods are merged

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
