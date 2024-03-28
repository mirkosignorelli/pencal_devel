#' Extract model fits from step 1 of PRC-LMM
#'
#' Summary function to extract the estimated fixed effect parameters and
#' variances of the random effects from an object fitted using `fit_lmms`
#' 
#' @param object the output of `fit_lmms`
#' @param yname a character giving the name of the longitudinal
#' variable for which you want to extract information
#' @param what one of the following: `'betas'` for the estimates 
#' of the regression coefficients; `'tTable'` for the usual T table
#' produced by `nlme`; `'variances'` for the estimates of 
#' the variances (and covariances) of the random effects and of the
#' variance of the error term
#' @param ... additional arguments
#' 
#' @return A vector containing the estimated fixed-effect parameters if 
#' `what = 'betas'`, the usual T table produced by `nlme` if `what = 'tTable'`,
#' or the estimated variance-covariance matrix of the random
#' effects and the estimated variance of the error if `what = 'variances'`
#' 
#' @import nlme
#' @export
#' 
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
#' 
#' @seealso \code{\link{fit_lmms}}

summary.lmmfit = function(object, yname, what = 'betas', ...) {
  what = match.arg(what, choices = c('betas', 'tTable', 'variances'))
  ynames = object$call.info$y.names
  if (! yname %in% ynames) {
    mess = paste(yname, 'is not one of the longitudinal variables in', 
                 deparse(substitute(object)))
    stop(mess)
  }
  else {
    requireNamespace('nlme')
    pos = which(object$call.info$y.names == yname)
    lmm = object$lmm.fits.orig[[pos]]
    if (what == 'betas') fixef(lmm)
    else if (what == 'tTable') summary(lmm)$tTable
    else if (what == 'variances') VarCorr(lmm)
    else {
      warning('what argument should be one of: betas, tTable, variances')
    }
  }
}

