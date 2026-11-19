# Part D2
rec = "agh" #placeholder TODO

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
  cd_train <- cluster_data(train$attended, model.matrix(~ age, data = train), train$practice)
  cd_test <- cluster_data(test$attended, model.matrix(~ age, data = test), test$practice)
  
  
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
s_glm

# tests/testthat/test-cross_validation.R
library(testthat)

# small simulated dataset shaped like gp
make_gp <- function(n_prac = 20, n_per = 15, seed = 1) {
  set.seed(seed)
  practice <- factor(rep(seq_len(n_prac), each = n_per))
  age <- rnorm(n_prac * n_per)
  u <- rnorm(n_prac, 0, 1)[practice]
  attended <- rbinom(length(age), 1, plogis(-0.1 + 0.8 * age + u))
  data.frame(attended, age, practice)
}

# a fake fitter that records what it was given, so we can test the CV logic
# without fitting anything
spy_fitter <- function() {
  calls <- list()
  f <- function(train, test, ...) {
    calls[[length(calls) + 1]] <<- list(train = train, test = test, dots = list(...))
    -nrow(test)                      # some deterministic negative score
  }
  list(f = f, calls = function() calls)
}

gp_test <- make_gp()

test_that("glm fitter returns a single finite, non-positive log-lik", {
  idx <- gp_test$practice %in% 1:15
  s <- fitter(gp_test[idx, ], gp_test[!idx, ], method = "glm")
  expect_length(s, 1)
  expect_true(is.finite(s))
  expect_lte(s, 0)
})

test_that("fitter rejects unknown methods", {
  expect_error(fitter(gp_test, gp_test, method = "banana"), "method must be one of")
})

test_that("cross_validation returns named cv and se", {
  set.seed(1)
  out <- cross_validation(gp_test, fitter)
  expect_named(out, c("cv", "se"))
  expect_true(all(is.finite(out)))
  expect_lte(out[["cv"]], 0)
  expect_gte(out[["se"]], 0)
})

test_that("no practice appears in both train and test", {
  spy <- spy_fitter()
  set.seed(1)
  cross_validation(gp_test, spy$f, folds = 5)
  for (cl in spy$calls()) {
    expect_length(intersect(unique(cl$train$practice), unique(cl$test$practice)), 0)
  }
})

test_that("every row is tested exactly once across folds", {
  spy <- spy_fitter()
  set.seed(1)
  cross_validation(gp_test, spy$f, folds = 5)
  calls <- spy$calls()
  expect_length(calls, 5)
  expect_equal(sum(sapply(calls, function(cl) nrow(cl$test))), nrow(gp_test))
})

test_that("extra arguments are forwarded to the fitter", {
  spy <- spy_fitter()
  set.seed(1)
  cross_validation(gp_test, spy$f, method = "agh", K = 7)
  dots <- spy$calls()[[1]]$dots
  expect_equal(dots$method, "agh")
  expect_equal(dots$K, 7)
})

test_that("same seed gives identical results", {
  set.seed(2026); a <- cross_validation(gp_test, fitter)
  set.seed(2026); b <- cross_validation(gp_test, fitter)
  expect_identical(a, b)
})

test_that("invalid folds are rejected", {
  expect_error(cross_validation(gp_test, fitter, folds = 1))
  expect_error(cross_validation(gp_test, fitter, folds = "5"))
})

#this runs and all tests pass as required just need s_rec
#set.seed(2026); s_rec <- cross_validation(gp, fitter, method = rec, K = 7)


