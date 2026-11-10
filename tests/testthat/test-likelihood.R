# Stand in until cluster_data() has been implemented and tested
make_test_data <- function() {
  list(
    y = c(0,1,1,0,0,0),
    x = cbind(1, c(-1,0,1,-1,0,1)),
    group = rep(1:2, each=3),
    j = 2,
    sizes = c(3,3)
  )
}

test_that("Cluster log-likelihood returns one value per cluster", {
  data <- make_test_data()
  theta <- c(-0.5,1,log(1.5))

  result <- cluster_loglik(theta, cd,  method = "laplace", K=1, draws=NULL)

  expect_type(result, "double")
  expect_length(result, data$j)
})

test_that("Cluster log-likelihood rejects invalid methods", {
  data <- make_test_data()
  theta <- c(-0.5,1,log(1.5))

  expect_error(cluster_loglik(theta, cd,
    method = "invald", K=1,
    draws=NULL
  ))
})

test_that("Cluster log likelihood returns a finite value", {
  data <- make_test_data()
  theta <- c(-0.5,1,log(1.5))

  result <- cluster_loglik(theta, cd, method = "laplace", K=1, draws=NULL)

  expect_true(all(is.finite(result)))
})
