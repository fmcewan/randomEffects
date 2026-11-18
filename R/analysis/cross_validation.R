# Part D2
fitter<- function(train, test,f="glm",K=NULL,draws=NULL) { #default is glm fit
  
  if (f=="glm"){
    fit <- glm(attended ~ age, family = binomial, data = train) #fit model
  }
  
  else {
    #check is one of required fns
    if (!f %in% (c("agh","gh","is","laplace"))){
      stop("f must be one of 'agh', 'gh', 'is', 'laplace'")
    }
    
    #This sits in my group analsis file
    fit <- fit_method(train,f,K,draws)
    
  }
    
  p <- predict(fit, test, type = "response") #predict model on test dataset
  sum(dbinom(test$attended, 1, p, log = TRUE)) #let log-lik be the error measure
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
s_glm

s_rec <- cross_validation(data=gp,fitter=fitter,f="agh",K=10)      #TODO

