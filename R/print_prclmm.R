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
#' @references Signorelli, M., Spitali, P., Tsonaka, R. (2021). Poisson-Tweedie 
#' mixed-effects model: a flexible approach for the analysis of longitudinal RNA-seq
#' data. Statistical Modelling, 21 (6), 520-545. URL: https://doi.org/10.1177/1471082X20936017
#' @seealso \code{\link{fit_prclmm}}, \code{\link{summary.prclmm}}
 
print.prclmm = function(x, digits = 4, ...) {
  temp = getinfo_step3(x)
  mod = temp$model_info
  dat = temp$data_info
  paste('Fitted model:', mod$fitted_model) |> cat(); cat('\n')
  paste('Penalty function used:', mod$penalty) |> cat(); cat('\n')
  
  paste('Sample size:', dat$n) |> cat(); cat('\n')
  paste('Number of events:', dat$n_ev) |> cat(); cat('\n')

  paste('Bootstrap optimism correction:', mod$cboc) |> cat(); cat('\n')
  
  show = round(temp$coefficients, digits)
  paste('Penalized likelihood estimates (rounded to ', digits, ' digits):', 
        sep = '') |> cat(); cat('\n')
  print(show) # cat( ) would remove the variable names!
}
