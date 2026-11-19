#'Gauss-Hermite function
#'
#'@description produces Gauss–Hermite nodes and weights for weight function e−x2, by the Golub–Welsch method.
#'Function based off the one from the Monte-Carlo integration lecture.
#'
#'@param n number of nodes and weights to produce
#'
#'@return a list containing the nodes and the weights
#'
#'@examples
#'gh <- gauss_hermite(5)
#'gh$nodes
#'
#'@export
gauss_hermite <- function(n) {             
  if (length(n) != 1 || !is.finite(n) || n < 1 || n != as.integer(n)){
    stop("n must be a valid positive integer")
  }

  k <- seq_len(n-1)
  #building matrix of 0s
  J <- matrix(0, n, n)
  #making J a symmetric tridiagonal Jacobi matrix
  J[cbind(k, k + 1)] <- J[cbind(k + 1, k)] <- sqrt(k / 2)
  
  e <- eigen(J, symmetric = TRUE)
  
  list(
    nodes = e$values, 
    weights = sqrt(pi) * e$vectors[1, ]^2)
}
