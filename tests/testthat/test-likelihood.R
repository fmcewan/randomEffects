# Example test data 
y <- c(0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0)
X <- cbind(intercept = 1, age = rep(c(-1, 0, 0.5, 1), 3))
group <- rep(1:3, each = 4)
data <- list(y = y, X = X, group = group, J = 3, sizes = c(4L, 4L, 4L))
theta <- c(-0.5,1,log(1.5))

test_that("Invalid method types are rejected in cluster_loglik", {
  expect_error(cluster_loglik(theta, data, method=NA, K=1, draws=NULL), "Invalid method name.")
  expect_error(cluster_loglik(theta, data, method=123, K=1, draws=NULL), "Invalid method name.")
  expect_error(cluster_loglik(theta, data, method=c("gh", "agh"), K=1, draws=NULL), "Invalid method name.")
  expect_error(cluster_loglik(theta, data, method=NULL, K=1, draws=NULL), "Invalid method name.")
})

test_that("Unknown method is rejected in cluster_loglik", {
  expect_error(cluster_loglik(theta, data, method="unknown", K=1, draws=NULL), "Method not yet implemented.")
})
