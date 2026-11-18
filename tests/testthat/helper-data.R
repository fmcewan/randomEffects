# Predefined data to use in tests 
y <- c(0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0)
X <- cbind(intercept = 1, age = rep(c(-1, 0, 0.5, 1), 3))
group <- rep(1:3, each = 4)
cd <- list(y = y, X = X, group = group, J = 3, sizes = c(4L, 4L, 4L))
theta <- c(-0.5,1,log(1.5))

# Larger data set for fitting
set.seed(2026)
fit_group <- rep(1:30, each = 10)
fit_X <- cbind(intercept = 1, x = rnorm(300))
fit_u <- rnorm(30, mean = 0, sd = 1)
fit_y <- rbinom(300, 1, plogis(drop(fit_X %*% c(-0.5, 1)) + fit_u[fit_group]))
fit_cd <- list(y = fit_y, X = fit_X, group = fit_group, J = 30L,sizes = rep(10L, 30))
