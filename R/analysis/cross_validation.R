# Part D2
fitter<- function(train, test,method="glm",K=NULL,draws=NULL) { 
  
  #baseline GLM model
  
  if (method=="glm"){
    fit <- glm(attended ~ age, family = binomial, data = train) #fit model
    p <- predict(fit, test, type = "response") #predict model on test dataset
    return (sum(dbinom(test$attended, 1, p, log = TRUE))) #let log-lik be the error measure
  }
  

  #check is one of required fns
  if (!method %in% (c("agh","gh","is","laplace"))){
    stop("method must be one of 'agh', 'gh', 'is', 'laplace'")
    }

   
  #change data
  cd_train <- make_cd(train)
  cd_test <- make_cd(test) #TODO
  
  
  # IS needs its own draws matrix for each dataset (one row per practice)
  draws_train <- draws_test <- NULL
  if (method == "is") {
    draws_train <- is_draws(cd_train, R = draws)
    draws_test  <- is_draws(cd_test,  R = draws)
  }
  
  #and fit on training data
  fit_train <- fit_rint(cd_train,method,K,draws_train)
  
  # 2. score held-out practices at those estimates (no refitting)
  return (marginal_loglik(fit_train$par, cd_test, method, K = K, draws = draws_test))
}



#adapted from lecture 12
cross_validation <- function(data,fitter,...,folds=5) {
  
  
  stopifnot(is.numeric(folds), folds >= 2)
  
  practices <- unique(data$practice)
  fold_of <- sample(rep(1:folds, length.out = length(practices)))   # fold per practice
  names(fold_of) <- practices
  data$fold <- fold_of[as.character(data$practice)]             # each patient gets their practice's fold
  
  fold_score <- numeric(folds)
  
  for (k in 1:folds) {
    train <- data[data$fold != k, ]
    test  <- data[data$fold == k, ]
    fold_score[k] <- fitter(train, test,...)       # fitter fits on train, returns held-out log-lik on test
  }
  c(cv=mean(fold_score), se = stats::sd(fold_score) / sqrt(folds))
}

set.seed(2026)
s_glm<-cross_validation(data=gp,fitter=fitter)

set.seed(2026); s_rec <- cross_validation(gp, fitter, method = "agh", K = 7)


