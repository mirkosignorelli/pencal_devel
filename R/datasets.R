#' pbc2 dataset
#'
#' This list contains data from the Mayo Clinic primary biliary cirrhosis (PBC)
#' study (1974-1984). It comprises two datasets, one with the survival and baseline covariates
#' and the other with the longitudinal measurements. The datasets are a 
#' rearrangement of the `pbc2` dataframe from the `joineRML` package that makes
#' them more suitable for analysis within `pencal`
#' 
#' @docType data
#' @keywords datasets
#' @usage data(pbc2data)
#'
#' @format The list contains two data frames:
#' 
#' \itemize{
#' \item \code{baselineInfo} contains the subject indicator `id`, information about
#' the survival outcome (`time` and `event`) and the covariates `baselineAge`, `sex`
#' and `treatment`;
#' \item \code{longitudinalInfo} contains the subject `id` and the repeated measurement 
#' data: `age` is the age of the individual at each visit, `fuptime` the follow-up time
#' (time on study), and `serBilir`, `serChol`, `albumin`, `alkaline`, `SGOT`,
#' `platelets` and `prothrombin` contain the value of each covariate at the 
#' corresponding visit
#' }
#' 
#' @author Mirko Signorelli
#' 
#' @examples 
#' data(pbc2data)
#' head(pbc2data$baselineInfo)
#' head(pbc2data$longitudinalInfo)
"pbc2data"
#'
#' A fitted PRC LMM
#'
#' This list contains a fitted PRC LMM, where the CBOCP is
#' computed using 50 cluster bootstrap samples. It is
#' used to reduce the computing time in the example of
#' the function \code{performance_prc}. The simulated dataset 
#' on which the model was fitted was landmarked at t = 2.
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
#' Signorelli, M. (2023). pencal: an R Package for the Dynamic 
#' Prediction of Survival with Many Longitudinal Predictors.
#' arXiv preprint: arXiv:2309.15600
#' 
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
#' the function \code{survpred_prcmlpmm}. The simulated dataset 
#' on which the model was fitted was landmarked at t = 2.
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