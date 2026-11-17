#Tests for gauss-hermite function

test_that("gh produces correctly formatted output",{
  expect_equal(5, length(gauss_hermite(5)$nodes))
  expect_equal(5, length(gauss_hermite(5)$weights))
})

test_that("gh produces appropriate error message for invalid input",{
  expect_error(gauss_hermite(0), regexp = "n must be a valid positive integer")
  expect_error(gauss_hermite(c(1,2)), regexp = "n must be a valid positive integer")
  expect_error(gauss_hermite(Inf), regexp = "n must be a valid positive integer")
  expect_error(gauss_hermite("5"), regexp = "n must be a valid positive integer")
})

test_that("nodes are symmetric around 0",{
  gh<-gauss_hermite(3)
  expect_true(all.equal(gh$nodes, -rev(gh$nodes)))
})


test_that("all weights are positive",{
  expect_true(all(gauss_hermite(4)$weights > 0))
})

#Tests for gh
example_theta <- c(0, -0.5, log(0.8))

example_cd <- list(
  y = c(0,0,1,0,0,0),
  X = matrix(c(1,1,1,1,1,1,0.1, 0.36, -0.2, 0, -0.23, -0.05), nrow = 6),
  group = c(1,1,2,2,3,3),
  J = 3,
  cluster_size = c(2,2,2)
)

invalid_cd <- list(
  y = c(0,0,"1",0,0,0),
  X = matrix(c(1,1,1,1,1,"one",0.1, 0.36, -0.2, 0, -0.23, -0.05), nrow = 6),
  group = c(1,1,2,2,3,3),
  J = 3,
  cluster_size = c(2,2,2)
)


void_cd <- list(
  y = c(0,0,1,0,0,NA),
  X = matrix(c(1,1,NA,1,1,1,0.1, 0.36, -0.2, 0, -0.23, -0.05), nrow = 6),
  group = c(1,1,2,2,3,3),
  J = 3,
  cluster_size = c(2,2,2)
)

test_that("gh produces correctly formatted output",{
  expect_true(length(gh(example_theta, example_cd, 4)) == 3)
})

test_that("gh produces errors for invalid inputs",{
  expect_error(gh(c(0, "-0.5", log(0.8)), example_cd, 4), regexp = "theta must contain exactly 3 numeric parameters")
  expect_error(gh(c(-0.5, log(0.8)), example_cd, 4), regexp = "theta must contain exactly 3 numeric parameters")
  expect_error(gh(example_theta, example_cd, 0), regexp = "K must be a valid positive integer")
  expect_error(gh(example_theta, invalid_cd, 3), regexp = "All columns in data must be numeric")
  expect_error(gh(example_theta, void_cd, 3), regexp = "Data must not contain NAs")
})

