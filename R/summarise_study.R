#' Performance measures for a single scenario
#'
#' Internal helper for [summarise_study()].
#'
#' @param results Rows of simulation results for one scenario.
#' @return A one-row data frame of performance measures and MCSEs.
#' @noRd
measure <- function(results) {
  
  #drop non-converged rows
  conv_rate<-mean(results$converged %in% TRUE) #store converged rows first
  results <- results[results$converged %in% TRUE, ]
  
  #change based on column names
  est<-results$est
  true<-results$true
  lower <-results$lower
  upper <- results$upper
  nsim<-length(est)
  
  #compute typical measure values along with MCSE
  
  #bias
  bias<- mean(est - true)
  
  #standard error
  s1<- est - mean(est)
  emp_se <- sqrt((1/(nsim-1)) * sum(s1^2) )
  
  #mean squared error
  s2 <- est - true
  mse <- 1/nsim * sum(s2^2)
  
  #coverage
  covers <- lower <= true & true <= upper
  coverage<-mean(covers)
  
  #mcse
  bias_mcse<-emp_se/sqrt(nsim)
  emp_se_mcse<-emp_se/sqrt(2*(nsim-1))
  mse_mcse<-stats::sd((est-true)^2)/sqrt(nsim)
  coverage_mcse <- sqrt(coverage * (1 - coverage) / nsim)
  
  
  data.frame(nsim = nsim, conv_rate = conv_rate, bias = bias, bias_mcse = bias_mcse,
             emp_se=emp_se,emp_se_mcse=emp_se_mcse,mse=mse,mse_mcse=mse_mcse,coverage=coverage,
             coverage_mcse=coverage_mcse
             )
}


#' Summarise simulation study results
#'
#' Computes standard performance measures for a simulation study, each with
#' its Monte Carlo standard error (MCSE), following Morris et al. (2019).
#' Results are summarised separately for every combination of number of
#' clusters, random-effect SD, method and parameter. Non-converged fits are
#' excluded from all measures except the convergence rate.
#'
#' @param results A data frame with one row per replicate x method x parameter,
#'   containing the columns:
#'   \describe{
#'     \item{rep}{Replicate number.}
#'     \item{m}{Number of clusters.}
#'     \item{sigma}{True random-intercept standard deviation.}
#'     \item{method}{Name of the fitting method.}
#'     \item{param}{Name of the parameter being estimated.}
#'     \item{true}{True value of the parameter.}
#'     \item{est}{Estimated value.}
#'     \item{se}{Model-based standard error of the estimate.}
#'     \item{lower, upper}{Lower and upper bounds of the 95 confidence interval.}
#'     \item{converged}{Logical; did the fit converge?}
#'     \item{time}{Time taken to fit, in seconds.}
#'   }
#'
#' @return A data frame with one row per combination of \code{m},
#'   \code{sigma}, \code{method} and \code{param}, containing:
#'   \describe{
#'     \item{nsim}{Number of converged replicates used.}
#'     \item{conv_rate}{Proportion of replicates that converged.}
#'     \item{bias, bias_mcse}{Bias, mean(est - true), and its MCSE.}
#'     \item{emp_se, emp_se_mcse}{Empirical standard error, sd(est), and its MCSE.}
#'     \item{mse, mse_mcse}{Mean squared error and its MCSE.}
#'     \item{coverage, coverage_mcse}{Proportion of 95 CIs containing the
#'       true value, and its MCSE.}
#'   }
#'
#' @references Morris, T. P., White, I. R. and Crowther, M. J. (2019).
#'   Using simulation studies to evaluate statistical methods.
#'   \emph{Statistics in Medicine}, 38(11), 2074--2102.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' make <- function(meth) {
#'   est <- rnorm(n, mean = 0.5, sd = 1)
#'   data.frame(rep = 1:n, m = 20, sigma = 1, method = meth,
#'              param = "beta1", true = 0.5, est = est, se = 1,
#'              lower = est - 1.96, upper = est + 1.96,
#'              converged = TRUE, time = 0.01)
#' }
#' results <- rbind(make("glm"), make("rint"))
#' summarise_study(results)
#'
#' @export
summarise_study <- function(results) {
  
  # check required columns (agreed with simulation output)
  req_names <- c("rep", "m", "sigma", "method", "param", "true", "est",
                 "se", "lower", "upper", "converged", "time")
  missing <- setdiff(req_names, names(results))
  
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # split into scenarios
  groups <- split(results, results[c("m", "sigma", "method", "param")], drop = TRUE)
  
  # measures for each scenario, with labels
  rows <- lapply(groups, function(g) {
    labels <- g[1, c("m", "sigma", "method", "param")]
    cbind(labels, measure(g))
  })
  
  # stack into one data frame
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
  }
  
  
  
