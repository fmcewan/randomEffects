# Laplace and adaptive Gauss-Hermite quadrature

#' Cluster log-likelihoods using Laplace approximation
#'
#' @param theta Vector of regression coefficients, one per column of cd$X,
#'   followed by log(sigma), where sigma is the standard deviation of the random intercepts
#' @param cd List containing:
#'   binary responses y,
#'   model matrix X,
#'   cluster indices group (containing every integer 1-J)
#'   and number of clusters J
#'
#' @details
#' Uses cluster_modes() for each cluster's mode and curvature, then applies the Laplace formula
#' (adapting laplace_j() from Practical 7, but using cluster_modes())
#' 
#' @return A named numeric vector of length J, 
#'   containing one approximate marginal log-likelihood per cluster, 
#'   ordered and labelled by cluster index. 
#'   Contributions are not summed across clusters.
#'
#' @examples
#' cd <- list(
#'  y = c(0, 1, 1, 0, 1, 0),
#'  X = cbind(intercept = rep(1, 6), x = c(-1, 0, 1, -1, 0, 1)),
#'  group = rep(1:2, each = 3), 
#'  J = 2
#' )
#' theta <- c(0.5, 0.8, log(1))
#' laplace_loglik(theta, cd)
#'  
laplace_loglik <- function(theta, cd) {
    
    # Storing/renaming values used later, for readability
    X <- cd$X
    y <- cd$y
    group <- cd$group 
    beta <- theta[seq_len(ncol(X))]
    sigma <- exp(theta[length(theta)])
    
    # Running cluster_modes and storing results
    modes <- cluster_modes(theta, cd)
    u_hat <- modes$u_hat
    H <- modes$H
    
    # Each patient's predictor (including practice effects)
    eta <- drop(X %*% beta) + u_hat[group]

    # Same logic to compute objective values as in cluster_modes
    softplus <- function(z) { pmax(z, 0) + log1p(exp(-abs(z)))} 
    log_prob <- ifelse(y == 1,-softplus(-eta),-softplus(eta))   
    h_hat <- drop(rowsum(log_prob,group,reorder=TRUE)) + dnorm(u_hat,mean=0,sd=sigma,log=TRUE)

    # From laplace_j (practical 7), slightly modified to avoid division
    log_lik <- h_hat + 0.5 * (log(2 * pi) - log(H))

    # Note for other group members: due to rowsum() this returned value will have labels attached
    # Use unname() in implementation if needed, I have kept it for consistency
    log_lik
}

#' Cluster log-likelihoods using adaptive Gauss-Hermite quadrature approximation
#' (adapted from Practical 7 extension)
#'
#' @param theta Vector of regression coefficients, 
#'   one per column of cd$X, 
#'   followed by log(sigma), where sigma is the standard deviation of the random intercepts
#' @param cd List containing:
#'   binary responses y, 
#'   model matrix X,
#'   cluster indices group and number of clusters J 
#' @param K Integer number of quadrature nodes per cluster, default=7
#' 
#' @details Uses cluster_modes() to centre and scale the nodes 
#'   supplied by gauss_hermite() for each cluster. 
#'   Weighted contributions are combined using log-sum-exp for numerical stability.
#'   (adapted from the adaptive Gauss-Hermite extension in Practical 7)
#'
#' @return A named vector of length cd$J, 
#'   containing one approximate marginal log-likelihood per cluster, 
#'   ordered and labelled by cluster index
#'
#' @examples
#' cd <- list(
#'  y = c(0, 1, 1, 0, 1, 0),
#'  X = cbind(intercept = rep(1, 6), x = c(-1, 0, 1, -1, 0, 1)),
#'  group = rep(1:2, each = 3), 
#'  J = 2
#' )
#' theta <- c(0.5, 0.8, log(1))
#' agh_loglik(theta, cd, K = 10)
# 
agh_loglik <- function(theta, cd, K = 7) {

    # Checking K is valid
    if (!is.numeric(K) || length(K) != 1 || !is.finite(K) || K < 1 || K != floor(K)) {
        stop("K must be a positive integer.")
    }

    # Storing/renaming for simplicity
    beta <- theta[seq_len(ncol(cd$X))]; sigma <- exp(theta[length(theta)])
    X <- cd$X; y <- cd$y
    group <- cd$group; J <- cd$J

    # Store modes and curvature of given parameters
    modes <- cluster_modes(theta, cd)
    u_hat <- unname(modes$u_hat); H <- unname(modes$H) 

    # Taking nodes and weights from gh (written by different group memeber, see gauss_hermite.R)
    gh <- gauss_hermite(K)
    nodes <- gh$nodes; weights <- gh$weights # Store in more convenient variables

    # Centering nodes at each mode and scale (adapted from practical 7 loglik_agq())
    scale <- 1 / sqrt(H) # Scale based on curvature
    # u is a JxK matrix (each row containing all nodes for that practice)
    u <- sweep(
        outer(sqrt(2) * scale, nodes),       # Multiplying nodes by the practice's scale
        MARGIN = 1, STATS = u_hat, FUN = "+" # Adding the mode to the nodes to recenter them
    )

    # Give every patient the integration points for their practice

    eta_fixed <- drop(X %*% beta)
    # eta has one row per patient and one column per node
    eta <- sweep(
        u[group, , drop = FALSE], # Selecting every patient's nodes (based on their practice)
        MARGIN = 1, STATS = eta_fixed, FUN = "+" # Adding that patient's fixed predictor to nodes
    )

    softplus <- function(z) { pmax(z, 0) + log1p(exp(-abs(z))) } # Same softplus

    # Repeating observed outcomes across node columns (one per patient)
    y_mat <- matrix(y, nrow = length(y), ncol = K)
    # Elementwise logprobability computation 
    log_prob <- ifelse(y_mat == 1, -softplus(-eta), -softplus(eta)) 

    # h contains one log-integrand value per practice and node
    h <- rowsum(log_prob, group, reorder = TRUE) + dnorm(u, mean = 0, sd = sigma, log = TRUE)

    # Adding corrected weights to log-integrand values (adapted from practical 7's loglik_agq())
    log_terms <- sweep(h, MARGIN = 2, STATS = log(weights) + nodes^2, FUN = "+")

    # Each row's maximum, without looping over practices (vector of length J)
    row_max <- log_terms[
        cbind(seq_len(J),  # generate J many (row, max_col) pairs into a Jx2 matrix
        max.col(log_terms, ties.method = "first")) # Each row's max column (no loop)
        ]

    # Rescaling based on the max entry for that practice (for numerical safety)
    shifted <- sweep(log_terms, MARGIN = 1, STATS = row_max, FUN = "-")
    log_sum <- row_max + log(rowSums(exp(shifted))) # Undoing rescaling

    # Finally, applying the change-of-variable factor (also adapted from practical 7)
    log_lik <- 0.5 * log(2) + log(scale) + log_sum

    setNames(log_lik, seq_len(J)) # Again, numbering results for consistency
}