#' Summary method for PRC-MLPMM model fits
#' 
#' @param object an object of class \code{prcmlpmm}
#' @param ... additional arguments
#' 
#' @return An object of class `sprcmlpmm`
#' @export
#' @author Mirko Signorelli
#' @references Signorelli, M., Spitali, P., Tsonaka, R. (2021). Poisson-Tweedie 
#' mixed-effects model: a flexible approach for the analysis of longitudinal RNA-seq
#' data. Statistical Modelling, 21 (6), 520-545. URL: https://doi.org/10.1177/1471082X20936017
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

