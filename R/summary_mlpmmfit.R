#' Extract model fits from step 1 of PRC-LMM
#'
#' Utility function to extract the MLPMM summaries from a model fit
#' obtained through `fit_mlpmms`
#' 
#' @param object the output of `fit_lmms`
#' @param yname a character giving the name of one of the longitudinal
#' outcomes modelled within one of the MLPMM
#' @param ... additional arguments
#' 
#' @return The model summary as returned by `summary.multlcmm`
#' 
#' @import lcmm
#' @export
#' 
#' @author Mirko Signorelli
#' @references 
#' Signorelli, M. (2024). pencal: an R Package for the Dynamic 
#' Prediction of Survival with Many Longitudinal Predictors. 
#' To appear in: The R Journal. Preprint: arXiv:2309.15600
#' 
#' Signorelli, M., Spitali, P., Al-Khalili Szigyarto, C, 
#' The MARK-MD Consortium, Tsonaka, R. (2021). 
#' Penalized regression calibration: a method for the prediction 
#' of survival outcomes using complex longitudinal and 
#' high-dimensional data. Statistics in Medicine, 40 (27), 6178-6196.
#' DOI: 10.1002/sim.9178
#' 
#' @seealso \code{\link{fit_mlpmms}} and \code{\link[lcmm]{summary.multlcmm}}

summary.mlpmmfit = function(object, yname, ...) {
  if (!is.character(yname)) stop('yname should be a character value')
  ynames = object$call.info$y.names
  check = yname %in% do.call(c, ynames)
  if (!check) {
    mess = paste(yname, 'is not one of the response variables in', 
                 deparse(substitute(object)))
    stop(mess)
  }
  else {
    requireNamespace('lcmm')
    pos = which( sapply(ynames, function(x) yname %in% x) )
    mlpmm = object$mlpmm.fits.orig[[pos]]
    summary(mlpmm)
  }
}

