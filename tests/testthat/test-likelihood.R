test_that("Cluster log-likelihood rejects invalid methods", {
  theta <- c(-0.5,1,log(1.5))

  expect_error(cluster_loglik(theta, cd, method = "invalid", K=1, draws=NULL))
})

test_that("Invalid method types are rejected in cluster_loglik", {

  expect_error(cluster_loglik(NA))
  expect_error(cluster_loglik(123))
  expect_error(cluster_loglik(c("gh", "agh")))
  expect_error(cluster_loglik(NULL))

})

