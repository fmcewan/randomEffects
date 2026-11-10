#------------Helper--------------#
#' Compute measures of data
#' 
#' @param
#' @return
#' @export
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
  mse_mcse<-sd((est-true)^2)/sqrt(nsim)
  coverage_mcse <- sqrt(coverage * (1 - coverage) / nsim)
  
  
  data.frame(nsim = nsim, conv_rate = conv_rate, bias = bias, bias_mcse = bias_mcse,
             emp_se=emp_se,emp_se_mcse=emp_se_mcse,mse=mse,mse_mcse=mse_mcse,coverage=coverage,
             coverage_mcse=coverage_mcse
             )
}


#' Summarise simulation study results
#'
#' @param results A data frame with columns:replicate number,number of clusters,size of cluster,true parameters,estimation method,estimated intercept,slop,standard deviation,standard errors,convergence signature possible(time,error)
#' @return A data frame with one row per scenario × method × param and returned estimates
#' @export
summarise_study <- function(results) {
  
  #Check column names are as expected from simulation study - this can be changed
  req_names<-c("replicate number","number of clusters","size of clusters","true parameter values","estimation method","estimated intercept","slope","standard deviation","standard errors","converged","time","error") #hard code
  missing<-setdiff(req_names,names(results))
  
  if (length(missing)>0){
    stop(paste(missing, collapse = ", "))
  }
  
  
  
}