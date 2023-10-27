#' Extract model fits from step 1 of PRC-LMM
#'
#' Utility function to extract the estimated fixed effect parameters and
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
#' 
#' @examples
#' # generate example data
#' set.seed(1234)
#' p = 4 # number of longitudinal predictors
#' simdata = simulate_prclmm_data(n = 100, p = p, p.relev = 2, 
#'              seed = 123, t.values = c(0, 0.2, 0.5, 1, 1.5, 2))
#'  
#' # specify options for cluster bootstrap optimism correction
#' # procedure and for parallel computing 
#' do.bootstrap = FALSE
#' # IMPORTANT: set do.bootstrap = TRUE to compute the optimism correction!
#' n.boots = ifelse(do.bootstrap, 100, 0)
#' more.cores = FALSE
#' # IMPORTANT: set more.cores = TRUE to speed computations up!
#' if (!more.cores) n.cores = 2
#' if (more.cores) {
#'    # identify number of available cores on your machine
#'    n.cores = parallel::detectCores()
#'    if (is.na(n.cores)) n.cores = 2
#' }
#' 
#' # step 1 of PRC-LMM: estimate the LMMs
#' y.names = paste('marker', 1:p, sep = '')
#' step1 = fit_lmms(y.names = y.names, 
#'                  fixefs = ~ age, ranefs = ~ age | id, 
#'                  long.data = simdata$long.data, 
#'                  surv.data = simdata$surv.data,
#'                  t.from.base = t.from.base,
#'                  n.boots = n.boots, n.cores = n.cores)
#'                  
#' # estimated betas and variances for the 5th marker:
#' getlmm(step1, 'marker2', 'betas')
#' getlmm(step1, 'marker2', 'tTable')
#' getlmm(step1, 'marker2', 'variances')

getlmm = function(object, yname, what = 'betas') {
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

