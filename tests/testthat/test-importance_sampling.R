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

# Test was AI suggested
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

test_that("is estimates lie within 4 of their own MC SE's of integrate()", {
  set.seed(10)
  cd <- sim_cd(20, 8, c(-0.5, 1), 1.5)
  theta <- c(-0.5, 1, log(1.5))
  est <- importance_sampling(theta, cd, is_draws(cd, R = 2000))
  exact <- exact_loglik(theta, cd)
  z <- (est - exact) / attr(est, "mc_se")     # error measured in SEs
  expect_true(all(abs(z) < 4))
  expect_lt(abs(mean(z)), 1)                  # no systematic offset
})

test_that("the same draws give exactly the same answer", {
  set.seed(11)
  cd <- sim_cd(15, 5, c(-0.5, 1), 1)
  Z <- is_draws(cd, R = 100)
  theta <- c(-0.5, 1, 0)
  expect_identical(importance_sampling(theta, cd, Z), importance_sampling(theta, cd, Z))
})

test_that("with fixed draws the estimate is smooth in theta", {
  set.seed(12)
  cd <- sim_cd(30, 5, c(-0.5, 1), 1)
  Z <- is_draws(cd, R = 100)
  theta <- c(-0.5, 1, 0)
  change <- sum(importance_sampling(theta + c(0, 0, 1e-5), cd, Z)) - sum(importance_sampling(theta, cd, Z))
  expect_lt(abs(change), 1e-3)  # a tiny step in theta gives a tiny change
})

# Test was AI suggested
test_that("the MC SE halves when R is quadrupled", {
  set.seed(13)
  cd <- sim_cd(30, 8, c(-0.5, 1), 1.5)
  theta <- c(-0.5, 1, log(1.5))
  se1 <- attr(importance_sampling(theta, cd, is_draws(cd, R = 1000)), "mc_se")
  se4 <- attr(importance_sampling(theta, cd, is_draws(cd, R = 4000)), "mc_se")
  expect_equal(median(se4 / se1), 0.5, tolerance = 0.15)
})

test_that("all-0 and all-1 clusters, and large sigma, give finite answers", {
  cd <- list(y = c(0,0,0,0, 1,1,1,1, 0,1,0,1), X = cbind(1, rep(c(-1, 0, 1, 2), 3)),
             group = rep(1:3, each = 4), J = 3)
  set.seed(14)
  est <- importance_sampling(c(0, 1, log(20)), cd, is_draws(cd, R = 200))
  expect_true(all(is.finite(est)))
  expect_true(all(is.finite(attr(est, "mc_se"))))
})

test_that("importance_sampling rejects missing or wrong-sized draws", {
  cd <- list(y = c(0, 1, 1, 0), X = cbind(1, c(-1, 0, 1, 2)), group = c(1, 1, 2, 2), J = 2)
  expect_error(importance_sampling(c(0, 1, 0), cd, NULL), "draws")
  expect_error(importance_sampling(c(0, 1, 0), cd, matrix(0, 5, 10)), "draws")
})