#Tests for cluster_data() function
#Example from roxygen comments
y_example <- c(0,0,0,1,1,0)
X_example <- matrix(c(1,1,1,1,1,1,0.1, 0.36, -0.2, 0, -0.23, -0.05), nrow = 6)
group_example <- factor(c(1,2,3,1,2,3))

test_that("cluster_data() works for valid inputs",{
  expect_true(all(cluster_data(y_example, X_example, group_example)$y == y_example))
  expect_true(all(cluster_data(y_example, X_example, group_example)$X == X_example))
  expect_true(all(cluster_data(y_example, X_example, group_example)$group == c(1,2,3,1,2,3)))
  expect_true(cluster_data(y_example, X_example, group_example)$J == 3)
  expect_true(all(cluster_data(y_example, X_example, group_example)$cluster_size == c(2,2,2)))
})

test_that("cluster_data() throws appropriate error when y is invalid",{
  expect_error(cluster_data(c(0,NA, 1, 1, 1, NA), X_example, group_example), regexp = "y must contain only numeric 0s and 1s")
  expect_error(cluster_data(c(2,3,2,3,3,3), X_example, group_example), regexp = "y must contain only numeric 0s and 1s")
})

test_that("cluster_data() throws appropriate error when X is invalid",{
  expect_error(cluster_data(y_example, matrix(c("one", "one", "one", "one", "one", "one",0.1, 0.36, -0.2, 0, -0.23, -0.05), nrow = 6), group_example), regexp = "X must be a numeric matrix with no NAs" )
  expect_error(cluster_data(y_example, c(c(1,1,1,1,1,1),c(0.1, 0.36, -0.2, 0, -0.23, -0.05)), group_example), regexp = "X must be a numeric matrix with no NAs" )
  expect_error(cluster_data(y_example, matrix(c(NA, NA, NA, NA, NA, NA,0.1, 0.36, -0.2, 0, -0.23, -0.05), nrow = 6), group_example), regexp = "X must be a numeric matrix with no NAs" )
  
})

test_that("cluster_data() throws appropriate error when y, X and group do not have an equal number of observations",{
  expect_error(cluster_data(c(0,1), X_example, group_example), regexp = "y, X and group must have the same number of observations")
  expect_error(cluster_data(y_example, matrix(c(1,1,1,1,1,0.1, 0.36, -0.2, 0, -0.23), nrow = 5), group_example), regexp = "y, X and group must have the same number of observations")
  expect_error(cluster_data(y_example, X_example, c(1,2,3,1)), regexp = "y, X and group must have the same number of observations")
})

test_that("cluster_data() throws appropriate error when group is not a factor or contains NAs",{
  expect_error(cluster_data(y_example, X_example, c(1,2,3,1,2,3)), regexp = "group must be a factor containing no NA values")
  expect_error(cluster_data(y_example, X_example, factor(c(1, 2, NA, NA, 2, 3))), regexp = "group must be a factor containing no NA values")
})

#Tests for rint_sim() function

J_test <- 3
m_test <- 4
beta_test <- c(1,0.2)
log_sigma_test <- log(0.8)
valid_test <- rint_sim(J_test, m_test, beta_test, log_sigma_test)

test_that("rint_sim() works for valid inputs",{
  expect_true(length(valid_test$y) == J_test*m_test)
  expect_true(all(dim(valid_test$X) == c(J_test*m_test, 2)))
  expect_true(length(valid_test$group) == J_test*m_test)
  expect_true(valid_test$J == J_test)
})

test_that("rint_sim() fails for invalid J",{
  expect_error(rint_sim("3", m_test, beta_test, log_sigma_test), regexp = "J must be an integer")
  expect_error(rint_sim(NA, m_test, beta_test, log_sigma_test), regexp = "J must be an integer")
  expect_error(rint_sim(3.2, m_test, beta_test, log_sigma_test), regexp = "J must be an integer")
})

test_that("rint_sim() fails for invalid m",{
  expect_error(rint_sim(J_test, "4", beta_test, log_sigma_test), regexp = "m must be an integer")
  expect_error(rint_sim(J_test, NA, beta_test, log_sigma_test), regexp = "m must be an integer")
  expect_error(rint_sim(J_test, 4.6, beta_test, log_sigma_test), regexp = "m must be an integer")
})

test_that("rint_sim() fails for invalid beta",{
  expect_error(rint_sim(J_test, m_test, c(1), log_sigma_test), regexp = "beta must contain exactly 2 numeric parameters")
  expect_error(rint_sim(J_test, m_test, c(NA, 3), log_sigma_test), regexp = "beta must contain exactly 2 numeric parameters")
  expect_error(rint_sim(J_test, m_test, c("1", "2"), log_sigma_test), regexp = "beta must contain exactly 2 numeric parameters")
})

test_that("rint_sim() fails for invalid log_sigma",{
  expect_error(rint_sim(J_test, m_test, beta_test, "2"), regexp = "log_sigma must be a numeric value")
})
