# Laplace and adaptive Gauss-Hermite quadrature

# laplace log likelihood, adapting laplace_j() from Practical 7 
# but using cluster_modes() instead of optimize()
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