#'Cluster Data function
#'
#'@description Prepares and validates the data before calculating the marginal log-likelihood
#'
#'@param y Vector containing 0s and 1s representing the binary response variable
#'@param X Matrix with one row per observation containing the covariates
#'@param group Factor which identifies which cluster each observation belongs to
#'
#'@return A list containing:
#'  - y : the response
#'  - X : the model matrix
#'  - group : a vector containing all distinct cluster labels
#'  - J : an integer corresponding to the number of clusters
#'  - cluster_size : the size of each cluster
#'  
#'@examples
#'y_example <- c(0,0,0,1,1,0)
#'X_example <- matrix(c(1,1,1,1,1,1,0.1, 0.36, -0.2, 0, -0.23, -0.05), nrow = 6)
#'group_example <- factor(c(1,2,3,1,2,3))
#'cluster_data(y_example, X_example, group_example)
#'
#'@export
cluster_data <- function(y, X, group){
  
  #Validate inputs
  if(anyNA(y) || !all(y == 0 | y == 1)){
    stop("y must contain only numeric 0s and 1s")
  }
  
  if (!is.matrix(X) || !is.numeric(X) || anyNA(X)){
    stop("X must be a numeric matrix with no NAs")
  }
  
  if(length(y) != nrow(X) || length(y) != length(group)){
    stop("y, X and group must have the same number of observations")
  }
  
  if (!is.factor(group) || anyNA(group)) {
    stop("group must be a factor containing no NA values")
  }
  
  #A vector containing all distinct cluster labels
  cluster_labels <- as.integer(droplevels(group))
  
  #An integer corresponding to the number of distinct cluster labels
  n_clusters <- length(unique(cluster_labels))
  
  #A vector stating the size of each cluster
  cluster_size <- tabulate(cluster_labels, nbins = n_clusters)
  
  return(list(
    y = y,
    X = X,
    group = cluster_labels,
    J = n_clusters,
    cluster_size = cluster_size))
}

#'rint_sim Function
#'
#'@description Simulate data according to the model with one standard-normal covariate
#'
#'@param J integer which represents the desired number of clusters
#'@param m integer which represents the desired number of observations within each cluster
#'@param beta vector of regression coefficients
#'@param log_sigma a float log(sigma) where sigma is the standard deviation of the random intercepts
#'
#'@return a list containing:
#'y = a vector of the simulated binary responses, 
#'X = a vector of the simulated design matrix, 
#'group = a vector assigning every observation to their cluster
#'J = an integer representing the number of clusters
#'
#'@examples
#'rint_sim(3, 4, c(1, 0.2), log(0.8))
#'
#'@export
rint_sim <- function(J, m, beta, log_sigma){
  #Checking J is a valid integer input
  #cannot check J is explicitly saved as an integer because often it saves it as numeric despite being a valid integer.
  if (!is.numeric(J) || J != floor(J)){
    stop("J must be an integer")
  }
  
  #Checking m is a valid integer input
  if (!is.numeric(m) || m != floor(m)){
    stop("m must be an integer")
  }

  if (length(beta) != 2 || !is.numeric(beta) || any(is.na(beta))) {
    stop("beta must contain exactly 2 numeric parameters")
  }

  if (!is.numeric(log_sigma)) {
    stop("log_sigma must be a numeric value")
  }
  
  #Generating the random intercepts for each cluster. 1 per cluster.
  u <- rnorm(J, mean = 0, sd = exp(log_sigma))
  
  #Generating the covariates. One covariate per observation.
  #mean = 0 and sd = 1 so utilising rnorm's default
  x <- rnorm(J*m)
  
  #Building the design matrix
  X <- cbind(intercept = 1, x = x)
  
  #Making a vector to later assign every observation to their appropriate random intercept.
  group <- rep(1:J, each = m) 
  
  #Calculating the probability of each observation having a response 1
  p <- plogis(drop(X %*% beta) + u[group])
  
  #Generating the binary reponses
  y <- rbinom(J*m, size = 1, prob = p)
  
  #returning an output 'cd' in line with expected input of other functions
  
  return(list(
    y = y, 
    X = X,
    group = group,
    J = J
  ))
}
