#' Compute the predicted survival probabilities obtained
#' from the PRC models
#'
#' This function computes the predicted survival probabilities 
#' for the for the PRC-MLPMM(U) and PRC-MLPMM(U+B) models proposed 
#' in Signorelli et al. (2021)
#' 
#' @param step2 the output of \code{\link{summarize_mlpmms}} 
#' (step 2 of the estimation of PRC-MLPMM)
#' @param step3 the output of \code{\link{fit_prcmlpmm}} (step 3 
#' of the estimation of PRC-MLPMM)
#' @param times numeric vector with the time points at which
#' to estimate the time-dependent AUC
#' 
#' @return A data frame with the predicted survival probabilities
#' computed at the supplied time points
#' 
#' @import foreach doParallel glmnet survival survivalROC survcomp
#' @export
#' 
#' @author Mirko Signorelli
#' @references 
#' Signorelli, M. (2024). pencal: an R Package for the Dynamic 
#' Prediction of Survival with Many Longitudinal Predictors. 
#' To appear in: The R Journal. Preprint: arXiv:2309.15600
#' 
#' Signorelli, M., Spitali, P., Al-Khalili Szigyarto, C, 
#' The MARK-MD Consortium, Tsonaka, R. (2021). 
#' Penalized regression calibration: a method for the prediction 
#' of survival outcomes using complex longitudinal and 
#' high-dimensional data. Statistics in Medicine, 40 (27), 6178-6196.
#' DOI: 10.1002/sim.9178
#' 
#' @seealso \code{\link{fit_mlpmms}} (step 1),
#' \code{\link{summarize_mlpmms}} (step 2) and 
#' \code{\link{fit_prcmlpmm}} (step 3).
#' 
#' @examples
#' data(fitted_prcmlpmm)
#'                    
#' # predict survival probabilities at times 3 to 6
#' surv.probs = survpred_prcmlpmm(fitted_prcmlpmm$step2, 
#'                  fitted_prcmlpmm$step3, times = 3:6)
#' ls(surv.probs)
#' head(surv.probs$predicted_survival)

survpred_prcmlpmm = function(step2, step3, times = 1) {
  call = match.call()
  # load namespaces
  requireNamespace('survival')
  requireNamespace('glmnet')
  
  ############################
  ##### CHECK THE INPUTS #####
  ############################
  if (!is.numeric(times)) stop('times should be numeric!')
  n.times = length(times)
  # checks on step 2 input
  temp = c('call', 'ranef.orig', 'n.boots')
  check1 = temp %in% ls(step2)
  mess1 = paste('step2 input should cointain:', do.call(paste, as.list(temp)) )
  if (sum(check1) != 3) stop(mess1)
  ranef.orig = step2$ranef.orig
  # checks on step 3 input
  temp = c('call', 'pcox.orig', 'surv.data', 'n.boots')
  check2 = temp %in% ls(step3)
  mess2 = paste('step2 input should cointain:', do.call(paste, as.list(temp)) )
  if (sum(check2) != 4) stop(mess2)
  baseline.covs = eval(step3$call$baseline.covs)
  pcox.orig = step3$pcox.orig
  surv.data = step3$surv.data
  n = length(unique(surv.data$id))
  
  ###############################
  ##### COMPUTE PREDICTIONS #####
  ###############################
  # reconstruct pieces  
  surv.orig = Surv(time = surv.data$time, event = surv.data$event)
  if (is.null(baseline.covs)) {
    X.orig = as.matrix(ranef.orig)
  }
  if (!is.null(baseline.covs)) {
    X0 = model.matrix(as.formula(baseline.covs), data = surv.data)
    X.orig = as.matrix(cbind(X0, ranef.orig))
    contains.int = '(Intercept)' %in% colnames(X.orig)
    if (contains.int) {
      X.orig = X.orig[ , -1] 
    }
  }
  beta.hat = coef(pcox.orig, s = 'lambda.min')
  temp = rownames(beta.hat)
  beta.hat = as.numeric(beta.hat)
  names(beta.hat) = temp
  # check if dimensionality correct
  if (ncol(X.orig) != length(beta.hat)) stop("Dimensions of beta.hat and X.train don't correspond")

  # convert glmnet Cox model to equivalent with survival package
  linpred.orig = X.orig %*% beta.hat
  f1 = as.formula(surv.orig ~ linpred.orig)
  cox.survival = coxph(f1, init = 1, 
     control = coxph.control(iter.max = 0))
  
  # compute survival probabilities
  temp.sfit = survfit(cox.survival,
       newdata = data.frame('linpred.orig' = linpred.orig),
                      se.fit = F, conf.int = F)
  shat.orig = t(summary(temp.sfit, 
                       times = times)$surv)
  preds = data.frame(rownames(ranef.orig),
                   shat.orig)
  names(preds) = c('id', paste('S(', times, ')', sep = '') )
  out = list('call' = call, 
             'predicted_survival' = preds)
  return(out)
}
