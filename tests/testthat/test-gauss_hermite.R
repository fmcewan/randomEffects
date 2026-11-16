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