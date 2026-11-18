# Part D1

#helper to fix data
make_cd <- function(data) {
  cluster_data(data$attended, model.matrix(~ age, data = data), droplevels(data$practice))
}

#set up
cd <- make_cd(gp)
set.seed(2026)
draws_is <- is_draws(cd,R=100)

#fit methods by fitrint
fits <- list(
  agh     = fit_rint(cd, "agh", K = 7),
  #gh      = fit_rint(cd, "gh",  K = 10),   # fails until gh_loglik exists #TODO
  is      = fit_rint(cd, "is",  draws = draws_is),
  laplace = fit_rint(cd, "laplace")
)

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



  
