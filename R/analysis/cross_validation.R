# Part D2
# Part D2
# # #Compare the random-intercept model with a pooled logistic regression,
# #glm(attended ~ age, family = binomial), 
# #by 5-fold cross-validation in which each fold holds out whole practices (Lecture 11). 
# Score each model by the log-likelihood of the held-out practices’ responses: 
#   for the random-intercept model, that is ∑jlogLj
# # over the held-out practices, at the parameters fitted to the other folds. 
# Report both cross-validated scores with their standard errors across folds. 
# Explain why holding out random patients instead would be wrong here [at most 100 words].

#data load
group <- 5   # replace with your group number
source("raw_data/group_data.R")
gp <- group_data(group)


#start by setting the glm model
glm_fit <- glm(attended~age,family=binomial,data=gp)

#we compare to random intercept model selected by sim study
##
##
##
##



#code adapted from lecture 12
set.seed(2026)

practices <- unique(gp$practice)                         # the 60 practices
fold_of <- sample(rep(1:K, length.out = length(practices)))  # 12 practices per fold



centre <- runif(J, 0, 10); u <- rnorm(J)          # group locations and effects
g <- rep(1:J, each = m)
xg <- centre[g] + rnorm(J * m, 0, 0.1)
yg <- 0.2 * xg + u[g] + rnorm(J * m, 0, 0.5)
nn1 <- function(xt, yt, xn) yt[sapply(xn, function(x0) which.min(abs(xt - x0)))]
cv_nn <- function(folds) mean(sapply(1:10, function(k) {
  tr <- folds != k
  mean((yg[!tr] - nn1(xg[tr], yg[tr], xg[!tr]))^2)
}))
random_folds <- sample(rep(1:10, length.out = J * m))
group_folds <- sample(rep(1:10, length.out = J))[g]
Jn <- 4000; cn <- runif(Jn, 0, 10)                 # new groups, to find the truth
xn <- cn + rnorm(Jn, 0, 0.1); yn <- 0.2 * xn + rnorm(Jn) + rnorm(Jn, 0, 0.5)
round(c(random_folds = cv_nn(random_folds), grouped_folds = cv_nn(group_folds),
        new_groups = mean((yn - nn1(xg, yg, xn))^2)), 2)