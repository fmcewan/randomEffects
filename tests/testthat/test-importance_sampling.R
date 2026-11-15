test_that("is_draws returns a J x R matrix with its df attached", {
  set.seed(1)
  Z <- is_draws(list(J = 5), R = 50, df = 4)
  expect_equal(dim(Z), c(5L, 50L))
  expect_equal(attr(Z, "df"), 4)
})

test_that("is_draws is reproducible with a seed", {
  set.seed(1); Z1 <- is_draws(list(J = 5), R = 50)
  set.seed(1); Z2 <- is_draws(list(J = 5), R = 50)
  Z3 <- is_draws(list(J = 5), R = 50) 
  expect_identical(Z1, Z2)
  expect_false(identical(Z1, Z3))
})

test_that("is_draws draws come from the t distribution", {
  set.seed(2)
  Z <- is_draws(list(J = 10), R = 10000, df = 5)
  # about 90% of t5 draws lie below its 90% quantile
  expect_equal(mean(Z < stats::qt(0.9, df = 5)), 0.9, tolerance = 0.01)
})

test_that("is_draws rejects invalid input", {
  cd <- list(J = 5)
  expect_error(is_draws(cd, R = 1), "at least 2")
  expect_error(is_draws(cd, R = 10.5), "whole number")
  expect_error(is_draws(cd, R = "a"), "whole number")
  expect_error(is_draws(cd, R = 10, df = 0), "positive")
  expect_error(is_draws(list(), R = 10), "J") 
})
