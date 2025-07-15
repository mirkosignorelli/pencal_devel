#' Print method for PRC MLPMM model fits
#' 
#' @param x an object of class \code{prcmlpmm} 
#' @param digits number of digits at which the printed estimated regression
#' coefficients should be rounded (default is 4)
#' @param ... additional arguments
#' 
#' @return Summary information about the fitted PRC-MLPMM model
#' 
#' @export
#' @author Mirko Signorelli
#' @references
#' Signorelli, M. (2024). pencal: an R Package for the Dynamic 
#' Prediction of Survival with Many Longitudinal Predictors. 
#' The R Journal, 16 (2), 134-153.
#' 
#' @seealso \code{\link{fit_prcmlpmm}}, \code{\link{summary.prcmlpmm}}
 
print.prcmlpmm = function(x, digits = 4, ...) {
  print.prclmm(x, digits)
}
