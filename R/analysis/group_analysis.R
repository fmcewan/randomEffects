# Part D1

#source data
library(here)
library(randomEffects)#everyone's functions
group <- 5
source(here("raw_data", "group_data.R"))
gp <- group_data(group)


#pull in four methods
methods <- c("agh","gh","is","laplace")
K <- c(7,10,1)
draws <- 100

recommended <- "agh" #TODO


#helper to fit method
fit_method <- function(data, method, K, draws = NULL) {
  
  X  <- model.matrix(~ age, data = data)          # intercept column + age
  cd <- cluster_data(data$attended, X, data$practice) #TODO
  
  negll <- function(theta) -marginal_loglik(theta, cd, method = method, K = K, draws = draws)
  
  time <- system.time(
    opt <- optim(c(0, 0, 0), negll, method = "BFGS")
  )[["elapsed"]]
  
  H <- numDeriv::hessian(negll, opt$par)
  V <- tryCatch(solve(H), error = function(e) matrix(NA, 3, 3))
  
  list(par = opt$par, se = sqrt(diag(V)),
       converged = opt$convergence == 0, time = time)
}

fits <- Map(function(m, k) fit_method(gp, m, k, draws = draws),
            methods, K)

names(fits) <- methods$method

#report estimates with CI

z <- qnorm(0.975)
results <- do.call(rbind,lapply(names(fits), function(m){
  
  #pull estimate out NB sigma is on log scale so need to sepearate
  est <- fits[[m]]$par 
  se <- fits[[m]]$se 
  lower  <- est - z*se
  upper <- est + z*se
  
  #build df
  return(data.frame(Method=m,
                    Parameter = c("beta0", "beta1", "sigma"),
                    est = c(est[1:2],exp(est[3])),
                    lower = c(lower[1:2],exp(lower[3])), 
                    upper = c(upper[1:2],exp(upper[3]))
  ))
}))

#interpret the effect of age as an odds ratio (do in write up): (this is beta 1 just exponentiate)



  
