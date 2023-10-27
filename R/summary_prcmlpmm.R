#' Summary method for PRC-MLPMM model fits
#' 
#' @param object an object of class \code{prcmlpmm}
#' @param ... additional arguments
#' 
#' @return An object of class `sprcmlpmm`
#' @export
#' @author Mirko Signorelli
#' @references 
#' Signorelli, M. (2023). pencal: an R Package for the Dynamic 
#' Prediction of Survival with Many Longitudinal Predictors.
#' arXiv preprint: arXiv:2309.15600
#' 
#' Signorelli, M., Spitali, P., Al-Khalili Szigyarto, C, 
#' The MARK-MD Consortium, Tsonaka, R. (2021). 
#' Penalized regression calibration: a method for the prediction 
#' of survival outcomes using complex longitudinal and 
#' high-dimensional data. Statistics in Medicine, 40 (27), 6178-6196.
#' DOI: 10.1002/sim.9178
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

