#' Visualize survival predictions for a fitted PRC model
#'
#' Visualize survival predictions for a fitted PRC model
#' 
#' @param step1 the output of \code{\link{fit_lmms}} or \code{\link{fit_mlpmms}}
#' @param step2 the output of \code{\link{summarize_lmms}} or
#' \code{\link{summarize_mlpmms}}
#' @param step3 the output of \code{\link{fit_prclmm}} or \code{\link{fit_prcmlpmm}}
#' 
#' @examples
#' # generate example data
#' simdata = simulate_prclmm_data(n = 100, p = 4, p.relev = 2, 
#'              t.values = c(0, 0.2, 0.5, 1, 1.5, 2),
#'              landmark = 2, seed = 123)
#'              
#' # estimate the PRC-LMM model
#' y.names = paste('marker', 1:p, sep = '')
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

survplot_prc = function(step1, step2, step3, ids, tmax = 5, res = 0.01,
                           legend.inset = -0.3, legend.space = 1) {
  par.init = par()
  
  x = seq(0, tmax, by = res)
  y = survpred_prclmm(step1, step2, step3, times = x)$predicted_survival
  y = subset(y, id %in% ids)
  
  par(mar = c(4, 4, 2, 7))
  plot(x, y[1, -1], ylim = c(0, 1), xlab = 't', ylab = 'S(t)',
       type = 'l', bty = 'L', col = 1)
  if (length(ids) > 1) {
    for (i in 2:nrow(y)) {
      points(x, y[i, -1], type = 'l', col = i)
    }
  }
  
  # add legend
  par(xpd = T)
  legend(x = "right", inset=c(legend.inset, 0), y$id, 
           title = 'subject', pch = 16,
           y.intersp = legend.space,
           col = 1:nrow(y), bty = 'n', cex = 1)
  # restore par values as they were
  par(bty = par.init$bty, xpd = par.init$xpd, mar = par.init$mar)
}

