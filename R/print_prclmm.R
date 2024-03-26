#' Print method for PRC-LMM model fits
#' 
#' @param x an object of class \code{prclmm} 
#' @param digits number of digits at which the printed estimated regression
#' coefficients should be rounded (default is 4)
#' @param ... additional arguments
#' 
#' @return Summary information about the fitted PRC-LMM model
#' 
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
#' @seealso \code{\link{fit_prclmm}}, \code{\link{summary.prclmm}}
 
print.prclmm = function(x, digits = 4, ...) {
  temp = getinfo_step3(x)
  mod = temp$model_info
  dat = temp$data_info
  paste('Fitted model:', mod$fitted_model) |> cat(); cat('\n')
  paste('Penalty function used:', mod$penalty) |> cat(); cat('\n')
  paste('Tuning parameters selected by CV:', temp$tuning) |> cat(); cat('\n')
  
  paste('Sample size:', dat$n) |> cat(); cat('\n')
  paste('Number of events:', dat$n_ev) |> cat(); cat('\n')

  paste('Bootstrap optimism correction:', mod$cboc) |> cat(); cat('\n')
  
  show = round(temp$coefficients, digits)
  paste('Penalized likelihood estimates (rounded to ', digits, ' digits):', 
        sep = '') |> cat(); cat('\n')
  print(show) # cat( ) would remove the variable names!
}
