#' Summary method for PRC MLPMM model fits
#' 
#' @param object an object of class \code{prcmlpmm}
#' @param ... additional arguments
#' 
#' @return An object of class `sprcmlpmm`
#' @export
#' @author Mirko Signorelli
#' @references 
#' Signorelli, M. (2024). pencal: an R Package for the Dynamic 
#' Prediction of Survival with Many Longitudinal Predictors. 
#' The R Journal, 16 (2), 134-153.
#' 
#' @seealso \code{\link{fit_prcmlpmm}}, \code{\link{print.prcmlpmm}}
 
summary.prcmlpmm = function(object, ...) {
  out = getinfo_step3(object)
  class(out) = 'sprcmlpmm'
  out
}

#' @method print sprcmlpmm   
#' @export

print.sprcmlpmm = function(x, digits = 4, ...) {
  print.sprclmm(x, digits)
}

