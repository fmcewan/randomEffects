# Part D1

library(here)
group <- 5
source(here("raw_data", "group_data.R"))
gp <- group_data(group)


#pull in four methods
methods <- data.frame(method = c("agh","gh","is","laplace"),
                      n = c(7,10,100,1))

recommended <- "agh" #TODO

#use marginal log-lik function #not made yet TODO



fit_method <- function(data,method,n){
  
  clusters <- split(data,data$practice)
  negll <- function(theta) -marginal_loglik(theta,clusters,method=method,n=n)
  
  time <- system.time(
    opt <- optim(c(0, 0, 0), negll, method = "BFGS")
  )[["elapsed"]]
  
  H <- numDeriv::hessian(negll, opt$par)
  
  V <- tryCatch(solve(H), error = function(e) matrix(NA, 3, 3))
  
  list(par = opt$par,                    # beta0, beta1, log sigma
       se = sqrt(diag(V)),
       converged = opt$convergence == 0,
       time = time)
}

fits <- lapply(seq_len(nrow(methods)), function(i)
  fit_method(gp, methods$method[i], methods$n[i]))
names(fits) <- methods$method
  
