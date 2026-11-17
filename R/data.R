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
#'  - n_clusters : an integer corresponding to the number of clusters
#'  - cluster_size : the size of each cluster
#'  
#'@examples
#'y_example <- c(0,0,0,1,1,0)
#'X_example <- matrix(c(1,1,1,1,1,1,0.1, 0.36, -0.2, 0, -0.23, -0.05), nrow = 6)
#'group_example <- c(1,2,3,1,2,3)
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
  
  if(!is.numeric(group) || anyNA(group)){
    stop("group must only contain integers")
  }
  
  #A vector containing all distinct cluster labels
  cluster_labels <- as.numeric(factor(group))
  
  #An integer corresponding to the number of distinct cluster labels
  n_clusters <- length(unique(cluster_labels))
  
  #A vector stating the size of each cluster
  cluster_size <- tabulate(group, nbins = n_clusters)
  
  return(list(
    y = y,
    X = X,
    group = cluster_labels,
    J = n_clusters,
    cluster_size = cluster_size))
}
