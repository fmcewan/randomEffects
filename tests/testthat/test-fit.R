peak <- c(-0.5, 1, 0)
fake_hill <- function(theta, ...) -sum((theta - peak)^2)
fake_flat <- function(theta, ...) -sum((theta[1:2] - peak[1:2])^2)

# TODO: Add integration tests once the real methods are merged 

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

  expect_equal(negative_loglik(theta, cd, "laplace", NULL, NULL),
               -fake_hill(theta))
})

