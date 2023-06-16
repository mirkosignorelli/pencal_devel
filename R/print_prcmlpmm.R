#' Print method for PRC-MLPMM model fits
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
#' @references Signorelli, M., Spitali, P., Tsonaka, R. (2021). Poisson-Tweedie 
#' mixed-effects model: a flexible approach for the analysis of longitudinal RNA-seq
#' data. Statistical Modelling, 21 (6), 520-545. URL: https://doi.org/10.1177/1471082X20936017
#' @seealso \code{\link{fit_prcmlpmm}}, \code{\link{summary.prcmlpmm}}
 
print.prcmlpmm = function(x, digits = 4, ...) {
  print.prclmm(x, digits)
}
