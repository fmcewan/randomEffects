# Part D1

library(here)
group <- 5
source(here("raw_data", "group_data.R"))
gp <- group_data(group)

#Functions needed

#data packages
source(here::here("R","data.R")) #cluster data
source(here::here("R","cluster_modes.R"))
source(here::here("R","fit.R")) #fitrint
source(here::here("R","likelihood.R")) #marginal log-lik

#functions packages
source(here::here("R","importance_sampling.R")) #is_draws + function
source(here::here("R","gauss_hermite.R"))
source(here::here("R","laplace_agh.R")) #agh

#helper to fix data
make_cd <- function(data) {
  cluster_data(data$attended, model.matrix(~ age, data = data), data$practice)
}

#set up
cd <- make_cd(gp)

#set seed and run random draws
set.seed(2026)
draws_is <- is_draws(cd,R=100)

#fit methods by fitrint
fits <- list(
  agh     = fit_rint(cd, "agh", K = 7),
  gh      = fit_rint(cd, "gh",  K = 10),  
  is      = fit_rint(cd, "is",  draws = draws_is),
  laplace = fit_rint(cd, "laplace")
)

#cluster modes still not converging

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
results


