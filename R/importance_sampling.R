#' Standardised t draws for importance sampling
#' 
#' Returns a J x R matrix of draws from a standard t distribution, one row per cluster. Later 
#' shifted and scaled to give importance sampling draws of each clusters random effect. The same
#' draws are used everytime the likelihood is evaluated
#' 
#' @param cd List of clustered data. Note only cd$J, the number of clusters is used
#' @param R Number of draws per cluster. Must be at least 2 as importance sampling estimate is
#' an average over draws, so need 2 or more to measure MC variability.
#' @param df Degrees of freedom of t distribution. Small df means heavy tails, so the proposal
#' covers the tails of the target safely & the importance weightings stay stable.
#' 
#' @return A `cd$J` by `R` numeric matrix of draws, with attribute `"df"`.
#' @export
#' 
#' @examples
#' cd <- list(J = 5)
#' set.seed(1)
#' Z <- is_draws(cd, R = 100, df = 5)
#' dim(Z)
#' attr(Z, "df")
is_draws <- function(cd, R, df = 5) {
  
  # Check the inputs
  if (!is.list(cd) || is.null(cd$J)) stop("`cd` must be a list containing `J` ")
  if (!is.numeric(R) || length(R) != 1 || !is.finite(R) || R < 2 || R != round(R))
    stop("`R` must be a single whole number of at least 2.")
  if (!is.numeric(df) || length(df) != 1 || is.na(df) || df <= 0)
    stop("`df` must be a single positive number.")
  
  # One row per cluster, R columns of standard t draws
  Z <- matrix(stats::rt(cd$J * R, df = df), nrow = cd$J, ncol = R)
  
  # Record df on the matrix, then return it
  attr(Z, "df") <- df
  Z
}