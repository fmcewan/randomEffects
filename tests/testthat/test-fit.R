peak <- c(-0.5, 1, 0)
fake_hill <- function(theta, ...) -sum((theta - peak)^2)
fake_flat <- function(theta, ...) -sum((theta[1:2] - peak[1:2])^2)

test_that("fit_rint recovers the peak from a known surface", {
  testthat::local_mocked_bindings(marginal_loglik = fake_hill)

  fit <- fit_rint(cd, method="laplace", start=c(0, 0, 0))
  
  expect_equal(fit$par, peak, tolerance=1e-3)
  expect_equal(fit$se, rep(sqrt(0.5), 3), tolerance=1e-4)
  expect_equal(fit$loglik, 0, tolerance=1e-6)
  expect_equal(fit$convergence, 0)
})

test_that("fit_rint returns NA errors on a flat surface", {
  testthat::local_mocked_bindings(marginal_loglik=fake_flat)

  fit <- fit_rint(cd, method="laplace", start=c(0,0,0))

  expect_true(all(is.na(fit$se)))
  expect_equal(fit$par[1:2], peak[1:2], tolerance=1e-3)
})

test_that("fit_rint uses logistic regression starting values without start", {
  testthat::local_mocked_bindings(marginal_loglik = fake_hill)

  fit <- fit_rint(cd, method = "laplace")

  expect_equal(fit$par, peak, tolerance = 1e-3)
  expect_equal(fit$convergence, 0)
})

test_that("negative_loglik is minus marginal_loglik", {
  testthat::local_mocked_bindings(marginal_loglik = fake_hill)

  expect_equal(negative_loglik(theta, cd, "laplace", NULL, NULL), -fake_hill(theta))
})

# Integration tests 
test_that("fit_rint with laplace converges on simulated data", {
  fit <- fit_rint(fit_cd, method = "laplace")

  expect_equal(fit$convergence, 0)
  expect_true(all(is.finite(fit$se)))
  expect_equal(fit$loglik, marginal_loglik(fit$par, fit_cd, "laplace"))
})

test_that("fit_rint with laplace agrees with lme4::glmer", {
  dat <- data.frame(y = fit_y, x = fit_X[, "x"], group = factor(fit_group))
  ref <- lme4::glmer(y ~ x + (1 | group), data = dat, family = binomial, nAGQ = 1)
  fit <- fit_rint(fit_cd, method = "laplace")

  expect_equal(unname(fit$par[1:2]), unname(lme4::fixef(ref)), tolerance = 1e-3)
  expect_equal(exp(fit$par[3]), lme4::getME(ref, "theta")[[1]], tolerance = 1e-3)
  expect_equal(fit$loglik, as.numeric(logLik(ref)), tolerance = 1e-4)
})

test_that("fit_rint with fixed IS draws converges on m = 20 data", {
  set.seed(4113)
  sim <- rint_sim(50, 20, c(-0.5, 1), log(1.5))
  cd <- cluster_data(sim$y, sim$X, factor(sim$group))
  fit <- fit_rint(cd, method = "is", draws = is_draws(cd, R = 100))
  expect_equal(fit$convergence, 0)
  expect_true(all(is.finite(fit$par)))
})
