#Example from roxygen comments
y_example <- c(0,0,0,1,1,0)
X_example <- matrix(c(1,1,1,1,1,1,0.1, 0.36, -0.2, 0, -0.23, -0.05), nrow = 6)
group_example <- c(1,2,3,1,2,3)

test_that("cluster_data works for valid inputs",{
  expect_true(all(cluster_data(y_example, X_example, group_example)$y == y_example))
  expect_true(all(cluster_data(y_example, X_example, group_example)$X == X_example))
  expect_true(all(cluster_data(y_example, X_example, group_example)$cluster_labels == c(1,2,3,1,2,3)))
  expect_true(cluster_data(y_example, X_example, group_example)$n_clusters == 3)
  expect_true(all(cluster_data(y_example, X_example, group_example)$cluster_size == c(2,2,2)))
})
