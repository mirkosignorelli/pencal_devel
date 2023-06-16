#' A fitted PRC LMM
#'
#' This list contains a fitted PRC LMM, where the CBOCP is
#' computed using 50 cluster bootstrap samples. It is
#' used to reduce the computing time in the example of
#' the function \code{performance_prc} 
#' 
#' @docType data
#' @keywords datasets
#' @usage data(fitted_prclmm)
#'
#' @format A list comprising step 2 and step 3 as obtained
#' during the estimation of a PRC LMM
#' 
#' @author Mirko Signorelli
#' @references 
#' Signorelli, M., Spitali, P., Al-Khalili Szigyarto, C, 
#' The MARK-MD Consortium, Tsonaka, R. (2021). 
#' Penalized regression calibration: a method for the prediction 
#' of survival outcomes using complex longitudinal and 
#' high-dimensional data. Statistics in Medicine.
#' DOI: 10.1002/sim.9178
#' 
#' @seealso \code{\link{performance_prc}}
#' 
#' @examples 
#' data(fitted_prclmm)
#' ls(fitted_prclmm)
"fitted_prclmm"
#'
#' A fitted PRC MLPMM
#'
#' This list contains a fitted PRC MLPMM. It is
#' used to reduce the computing time in the example of
#' the function \code{survpred_prcmlpmm} 
#' 
#' @docType data
#' @keywords datasets
#' @usage data(fitted_prclmm)
#'
#' @format A list comprising step 2 and step 3 as obtained
#' during the estimation of a PRC MLPMM
#' 
#' @author Mirko Signorelli
#' @references 
#' Signorelli, M., Spitali, P., Al-Khalili Szigyarto, C, 
#' The MARK-MD Consortium, Tsonaka, R. (2021). 
#' Penalized regression calibration: a method for the prediction 
#' of survival outcomes using complex longitudinal and 
#' high-dimensional data. Statistics in Medicine, 40 (27), 6178-6196.
#' DOI: 10.1002/sim.9178
#' 
#' @seealso \code{\link{survpred_prcmlpmm}}
#' 
#' @examples 
#' data(fitted_prcmlpmm)
#' ls(fitted_prcmlpmm)
"fitted_prcmlpmm"