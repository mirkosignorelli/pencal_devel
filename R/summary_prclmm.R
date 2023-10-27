#' Summary method for PRC-LMM model fits
#' 
#' @param object an object of class \code{prclmm} 
#' @param ... additional arguments
#' 
#' @return An object of class `sprclmm`
#' @importFrom methods is
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
#' @seealso \code{\link{fit_prclmm}}, \code{\link{print.prclmm}}
 
summary.prclmm = function(object, ...) {
  out = getinfo_step3(object)
  class(out) = 'sprclmm'
  out
}

#' @method print sprclmm    
#' @export

print.sprclmm = function(x, digits = 4, ...) {
  mod = x$model_info
  dat = x$data_info
  paste('Fitted model:', mod$fitted_model) |> cat(); cat('\n')
  paste('Penalty function used:', mod$penalty) |> cat(); cat('\n')
  
  paste('Sample size:', dat$n) |> cat(); cat('\n')
  paste('Number of events:', dat$n_ev) |> cat(); cat('\n')
  
  paste('Bootstrap optimism correction:', mod$cboc) |> cat(); cat('\n')
  
  show = round(x$coefficients, digits)
  paste('Penalized likelihood estimates (rounded to ', digits, ' digits):', 
        sep = '') |> cat(); cat('\n')
  print(show) # cat( ) would remove the variable names!
}

getinfo_step3 = function(x) {
  type = ifelse(is(x, 'prclmm'), 'PRC-LMM', 'PRCMLPMM')
  model_info = list(
    fitted_model = type,
    penalty = eval(x$call$penalty),
    cboc = ifelse(x$n.boots == 0, 'not computed', 
                  paste('computed using', x$n.boots, 'bootstrap samples'))
  )
  data_info = list(
    n = x$surv.data |> nrow(),
    n_events = x$surv.data$event |> sum()
  )
  coefficients = coef(x$pcox.orig, s = 'lambda.min') |> 
    as.matrix() |> t()
  out = list(model_info = model_info, 
             data_info = data_info, 
             coefficients = coefficients)
  out
}
