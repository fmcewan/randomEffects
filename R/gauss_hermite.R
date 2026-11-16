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
#'@param data The actual patient data used to evaluate the Bernoulli likelihood at each quadrature node
#'@param K number of nodes
#'
#'@examples
#'example_theta <- c(-0.5, log(0.8))
#'example_data <- data.frame(
#'   practice = c(1,1,2,2,3,3),
#'   age = c(-0.36, -2.03, 1.3, -0.89, -0.06, -0.05),
#'   attended = c(0,0,1,0,0,0)
#')
#'gh(example_theta, example_data, 4)
#'
#'
#'@export
gh <- function(theta, data, K){
  beta <- theta[1]
  sigma <- exp(theta[2])
  
  cluster <- factor(data$practice)
  n_clusters <- nlevels(cluster)
  
  x <- data$age
  y <- data$attended
  n <- length(y)
  
  gh <- gauss_hermite(K)
  
  #let u=sqrt(2)*sigma*x to transform the normal density into the required Gauss-Hermite weight function
  u <- sqrt(2)*sigma*gh$nodes
  
  #calculating beta*x +u for every quadrature node K
  #each row is an observation. beta*x will be consistent across the row.
  #each column corresponds to a quadrature node. u will be consistent down a column.
  linear_predictors <- outer(beta*x, rep(1,K)) + matrix(rep(u, each = n), nrow = n)
  
  #Calculate Bernoulli log probabilities. 
  #Note: 1-plogis(t) = plogis(-t)
  log_prob <- y*plogis(linear_predictors, log.p=TRUE)+ (1-y)*plogis(-linear_predictors, log.p=TRUE)
  
  #Sum log probabilities within each cluster (ie across all rows of the same cluster)
  #Repeats the calculation for each quadrature node
  cluster_loglikelihood <- sapply(seq_len(K), function(k){
    rowsum(log_prob[,k], cluster)[,1] #[,1] extracts the first (and only column) fo the result returned by rowsum()
  })
  
  #Take the logs of both weights and sqrt(pi) so operations can be performed on the log-scale
  #This reduces numerical underflow
  log_weights <- log(gh$weights/sqrt(pi))
  
  log_likelihoods <- apply(
    #adds the log weights to the log-likelihoods
    cluster_loglikelihood + matrix(log_weights, 
                                   nrow = n_clusters,
                                   ncol = K, 
                                   byrow = TRUE),
    #apply the function to each row of the matrix (ie every cluster)
    1, 
    #the function to apply is the log-sum-exp as suggested in the brief
    function(z){
      m <- max(z)
      m + log(sum(exp(z-m)))
    }
  )
  #return a vector containing the log likelihood for each cluster j
  return(log_likelihoods)
}
