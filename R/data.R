# cluster_data and rint_sim 

#'Cluster Data function
#'
#'@description Prepares and validates the data before calculating the marginal log-likelihood
#'
#'@param y Vector containing 0s and 1s representing the binary response variable
#'@param X Matrix with one row per observation containing the covariates
#'@param group Vector which identifies which cluster each observation belongs to
#'
#'@return A list containing:
#'  - y : the response
#'  - X : the model matrix
#'  - cluster_labels : a vector containing all distinct cluster labels
#'  - n_clusters : the number of clusters
#'  - cluster_size : the size of each cluster
#'  
#'@examples
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
  
  if(!is.integer(group) || anyNA(group)){
    stop("group must only contain integers")
  }
}