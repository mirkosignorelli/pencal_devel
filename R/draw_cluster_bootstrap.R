#' Draw a cluster bootstrap sample from a data frame in long format
#'
#' This function is part of the cluster bootstrap optimism correction
#' procedure described in Signorelli et al. (2021). Note 
#' that the function does not perform the random sampling, but it
#' extracts the correct records from a dataframe, given the ids of
#' the sampled clusters (subjects)
#' 
#' @param df a data frame in long format
#' @param idvar name of the subject id in \code{df} (it should be a
#' numeric id that ranges from 1 to n, without skipping values)
#' @param boot.ids identifiers of the subjects to be sampled
#' 
#' @return A data frame containing the bootstrapped observations
#' 
#' @import foreach
#' @keywords internal
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

draw_cluster_bootstrap = function(df, idvar, boot.ids) {
  # df = data frame
  # idvar = numeric id that ranges from 1 to n (no numbers skipped!)
  # boot.ids = ids of bootstrapped units
  requireNamespace('foreach')
  # avoid R CMD check note "no visible binding for global variable"
  i = NULL
  # proceed
  id.long = df[ , deparse(substitute(idvar))]
  n = max(id.long)
  if (length(boot.ids) != n) stop('bootstrap sample not of size n')
  is.first = !duplicated(boot.ids)
  # change the format of the new.ids variable
  # new.ids = ifelse(is.first, boot.ids, (n+1):(2*n))
  ratio = ifelse(n < 1000, 1000, 100000)
  if (n >= 100000) ratio = 1e8
  new.ids = ifelse(is.first, boot.ids, boot.ids + (1:n)/(ratio))
  # plot(boot.ids, new.ids)
  out = foreach(i = 1:n, .combine = 'rbind') %do% {
    keep = which(id.long == boot.ids[i])
    sub = df[keep, ]
    sub$new.id = new.ids[i]
    sub
  }
  return(out)
}
