# Predefined data to use in tests 
y <- c(0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0)
X <- cbind(intercept = 1, age = rep(c(-1, 0, 0.5, 1), 3))
group <- rep(1:3, each = 4)
cd <- list(y = y, X = X, group = group, J = 3, sizes = c(4L, 4L, 4L))
theta <- c(-0.5,1,log(1.5))
