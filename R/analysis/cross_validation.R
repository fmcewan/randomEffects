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


f_glm <- function(train, test) {
  fit <- glm(attended ~ age, family = binomial, data = train) #fit model
  p <- predict(fit, test, type = "response") #predict model on test dataset
  sum(dbinom(test$attended, 1, p, log = TRUE)) #let log-lik be the error measure
}

#we compare to random intercept model selected by sim study
##
##
##
##


#adapted from lecture 12

cross_validation <- function(data, K = 5,f) {
  
  practices <- unique(data$practice)
  fold_of <- sample(rep(1:K, length.out = length(practices)))   # fold per practice
  names(fold_of) <- practices
  data$fold <- fold_of[as.character(data$practice)]             # each patient gets their practice's fold
  
  fold_score <- numeric(K)
  for (k in 1:K) {
    train <- data[data$fold != k, ]
    test  <- data[data$fold == k, ]
    fold_score[k] <- f(train, test)       # f fits on train, returns held-out log-lik on test
  }
  c(cv=mean(fold_score), se = stats::sd(fold_score) / sqrt(K))
}

set.seed(2026)
s_glmcross_validation(data=gp,f=f_glm)

