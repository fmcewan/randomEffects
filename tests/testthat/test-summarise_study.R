make_fake <- function(meth = "glm", n = 2000) {
  est <- rnorm(n, mean = 2, sd = 1)
  data.frame(rep = 1:n, m = 20, sigma = 1, method = meth, param = "beta1",
             true = 2, est = est, se = 1, lower = est - 1.96, upper = est + 1.96,
             converged = TRUE, time = 0.1)
}

test_that("recovers known values on N(true, 1) data", {
  set.seed(2)
  out <- summarise_study(make_fake())
  expect_lt(abs(out$bias), 3 * out$bias_mcse)
  expect_lt(abs(out$emp_se - 1), 3 * out$emp_se_mcse)
  expect_lt(abs(out$coverage - 0.95), 3 * out$coverage_mcse)
})
