#' Summary for step 2 of PRC
#'
#' Summary function to extract basic descriptives from `summarize_lmms`
#' and `summarize_mlpmms`
#' 
#' @param object the output of `summarize_lmms` or `summarize_mlpmms`
#' @param ... additional arguments
#' 
#' @return Information about number of predicted random effects and sample size
#' 
#' @export
#' 
#' @author Mirko Signorelli
#' @references 
#' Signorelli, M. (2024). pencal: an R Package for the Dynamic 
#' Prediction of Survival with Many Longitudinal Predictors. 
#' The R Journal, 16 (2), 134-153.
#' 
#' @seealso \code{\link{summarize_lmms}}, \code{\link{summarize_mlpmms}}

summary.ranefs = function(object, ...) {
  dims = dim(object$ranef.orig)
  paste('Number of predicted random effect variables:', dims[2]) |> cat(); cat('\n')
  paste('Sample size:', dims[1]) |> cat()
}

