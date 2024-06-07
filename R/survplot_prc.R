#' Visualize survival predictions for a fitted PRC model
#'
#' Visualize survival predictions for a fitted PRC model
#' 
#' @param step1 the output of \code{\link{fit_lmms}} or \code{\link{fit_mlpmms}}
#' @param step2 the output of \code{\link{summarize_lmms}} or
#' \code{\link{summarize_mlpmms}}
#' @param step3 the output of \code{\link{fit_prclmm}} or \code{\link{fit_prcmlpmm}}
#' @param ids a vector with the identifiers of the subjects to show in the plot
#' @param tmax maximum prediction time to consider for the chart. Default is 5
#' @param res resolution at which to evaluate predictions for the chart. Default is 0.01
#' @param lwd line width
#' @param lty line type
#' @param legend.title legend title
#' @param legend.inset moves legend more to the left / right (default is -0.3)
#' @param legend.space interspace between lines in the legend (default is 1)
#' 
#' @importFrom graphics legend par points
#' @export
#' 
#' @author Mirko Signorelli
#' 
#' @references 
#' Signorelli, M. (2024). pencal: an R Package for the Dynamic 
#' Prediction of Survival with Many Longitudinal Predictors. 
#' To appear in: The R Journal. Preprint: arXiv:2309.15600
#' 
#' @examples
#' # generate example data
#' simdata = simulate_prclmm_data(n = 100, p = 4, p.relev = 2, 
#'              t.values = c(0, 0.2, 0.5, 1, 1.5, 2),
#'              landmark = 2, seed = 123)
#'              
#' # estimate the PRC-LMM model
#' y.names = paste('marker', 1:4, sep = '')
#' step1 = fit_lmms(y.names = y.names, 
#'                  fixefs = ~ age, ranefs = ~ age | id, 
#'                  long.data = simdata$long.data, 
#'                  surv.data = simdata$surv.data,
#'                  t.from.base = t.from.base,
#'                  n.boots = 0)
#' step2 = summarize_lmms(object = step1)
#' step3 = fit_prclmm(object = step2, surv.data = simdata$surv.data,
#'                    baseline.covs = ~ baseline.age,
#'                    penalty = 'ridge')
#'
#' # visualize the predicted survival for subjects 1, 3, 7 and 13                    
#' survplot_prc(step1, step2, step3, ids = c(1, 3, 7, 13), tmax = 6)

survplot_prc = function(step1, step2, step3, ids, tmax = 5, res = 0.01, lwd = 1, lty = 1,
                        legend.title = 'Subject', legend.inset = -0.3, legend.space = 1) {
  # fix for 'no visible binding for global variable ‘id’' note
  id = NULL
  par.init = par()
  
  x = seq(0, tmax, by = res)
  y = survpred_prclmm(step1, step2, step3, times = x)$predicted_survival
  y = subset(y, id %in% ids)
  
  par(mar = c(4, 4, 2, 7))
  plot(x, y[1, -1], ylim = c(0, 1), xlab = 't', ylab = 'S(t)',
       type = 'l', bty = 'L', col = 1, lwd = lwd, lty = lty)
  if (length(ids) > 1) {
    for (i in 2:nrow(y)) {
      points(x, y[i, -1], type = 'l', col = i, lwd = lwd, lty = lty)
    }
  }
  
  # add legend
  par(xpd = T)
  legend(x = "right", inset=c(legend.inset, 0), y$id, 
           title = legend.title, pch = 16,
           y.intersp = legend.space,
           col = 1:nrow(y), bty = 'n', cex = 1)
  # restore par values as they were
  par(bty = par.init$bty, xpd = par.init$xpd, mar = par.init$mar)
}

