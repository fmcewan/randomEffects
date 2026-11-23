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

#'Non-adaptive numerical integral based on Gauss-Hermite
#'
#'@description K-node Gauss-Hermite rule for the N(0,sigma^2) distribution of u
#'
#'@param theta Vector containing reparameterised beta and log sigma
#'@param cd List containing: 
#'   y: a numeric vector of binary responses,
#'   X: a finite numeric model matrix with one row per response,
#'   group: Holds integer identifiers (from 1 to J) matching patient to practice,
#'   J: integer stating the number of clusters 
#'   and cluster_size: the number of observations per cluster
#'@param K number of nodes
#'
#'@return a vector containing the log-likelihood for each cluster j
#'
#'@examples
#'example_theta <- c(0, -0.5, log(0.8))
#'example_cd <- list(
#' y = c(0,0,1,0,0,0),
#' X = matrix(c(1,1,1,1,1,1,0.1, 0.36, -0.2, 0, -0.23, -0.05), nrow = 6),
#' group = c(1,1,2,2,3,3),
#' J = 3,
#' cluster_size = c(2,2,2)
#')
#'example_K <- 4
#'gh(example_theta, example_cd, example_K)
#'
#'
#'@export
gh <- function(theta, cd, K){
  if (length(theta) != 3 || !is.numeric(theta)) {
    stop("theta must contain exactly 3 numeric parameters")
  }
  beta <- theta[1:2]
  sigma <- exp(theta[3])
  
  if (length(K) != 1 || !is.finite(K) || K < 1 || K != as.integer(K)){
    stop("K must be a valid positive integer")
  }
  
  if (!all(sapply(cd$y, is.numeric)) || !all(sapply(cd$X, is.numeric)) || !all(sapply(cd$group, is.numeric))) {
    stop("All columns in data must be numeric")
  }
  
  if (any(is.na(cd$y)) || any(is.na(cd$X)) ||any(is.na(cd$group))){
    stop("Data must not contain NAs")
  }
  
  cluster <- factor(cd$group)
  n_clusters <- cd$J
  
  x <- cd$X
  y <- cd$y
  n <- length(y)
  
  gh <- gauss_hermite(K)
  
  #let u=sqrt(2)*sigma*x to transform the normal density into the required Gauss-Hermite weight function
  u <- sqrt(2)*sigma*gh$nodes
  
  #calculating beta0*x0 + beta1*x1 +u for every quadrature node K 
  #In our case x0 will just be a column of 1s as we are working with one co-variate. 
  #The equation becomes beta0 + beta1*x1
  #each row is an observation. beta0 + beta1*x will be consistent across the row.
  #each column corresponds to a quadrature node. u will be consistent down a column.
  linear_predictors <- outer(as.vector(x %*% beta), rep(1,K)) + matrix(rep(u, each = n), nrow = n)
  
  #Calculate Bernoulli log probabilities. 
  #Note: 1-plogis(t) = plogis(-t)
  log_prob <- y*plogis(linear_predictors, log.p=TRUE)+ (1-y)*plogis(-linear_predictors, log.p=TRUE)
  
  #Sum log probabilities within each cluster (ie across all rows of the same cluster)
  #Repeats the calculation for each quadrature node
  cluster_loglikelihood <- rowsum(log_prob, cluster)
  
  #Take the logs of both weights and sqrt(pi) so operations can be performed on the log-scale
  #This reduces numerical underflow
  log_weights <- log(gh$weights/sqrt(pi))
  
  log_likelihoods <- row_logsumexp(
    #adds the log weights to the log-likelihoods
    cluster_loglikelihood + matrix(
      log_weights, 
      nrow = n_clusters,
      ncol = K, 
      byrow = TRUE)
    )
  #return a vector containing the log likelihood for each cluster j
  return(log_likelihoods)
}


#'Function to calculate the true marginal likelihood of a cluster
#'
#'Used for testing the accuracy of gh in test-gauss_hermite.R
#'
#'@param theta Vector containing reparameterised beta and log sigma
#'@param cd List containing: 
#'   y: a numeric vector of binary responses,
#'   X: a finite numeric model matrix with one row per response,
#'   group: Holds integer identifiers (from 1 to J) matching patient to practice,
#'   J: integer stating the number of clusters 
#'   and cluster_size: the number of observations per cluster
#'@param j cluster to integrate
#'
#'@return log marginal likelihood of the given cluster as a float
#'
cluster_integral <- function(theta, cd, j){
  beta <- theta[1:2]
  sigma <- exp(theta[3])
  
  indices <- cd$group == j
  
  y <- cd$y[indices]
  X <- cd$X[indices, , drop=FALSE] #extracts the rows relevant to the specific cluster. 
  #drop == FALSE keeps X as a matrix
  
  #Integrate over the random effect u
  result <- stats::integrate(
    function(u){
      linear_predictor <- outer(as.vector(X%*%beta),u, "+")
      #each row of linear_predictor is an observation and each column represents a different random effect
      
      #log_lik stores the log_likelihood for each value of u 
      log_lik <- colSums(
        y *plogis(linear_predictor, log.p = TRUE)+ (1-y)*plogis(-linear_predictor, log.p = TRUE)
      )
      
      exp(log_lik) * dnorm(u, mean = 0, sd = sigma)
    },
    lower = Inf,
    upper = Inf
  )
  log(result$value)
}
