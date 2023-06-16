#' Predictive performance of the PRC-LMM and PRC-MLPMM models
#'
#' This function computes the naive and optimism-corrected
#' measures of performance (C index and time-dependent AUC) 
#' for the PRC models proposed 
#' in Signorelli et al. (2021). The optimism
#' correction is computed based on a cluster bootstrap
#' optimism correction procedure (CBOCP)
#' 
#' @param step2 the output of either \code{\link{summarize_lmms}} 
#' or \code{\link{summarize_mlpmms}} (step 2 of the estimation of
#' PRC)
#' @param step3 the output of \code{\link{fit_prclmm}} or
#' \code{\link{fit_prcmlpmm}} (step 3 of PRC)
#' @param times numeric vector with the time points at which
#' to estimate the time-dependent AUC
#' @param n.cores number of cores to use to parallelize part of
#' the computations. If \code{ncores = 1} (default), 
#' no parallelization is done. Pro tip: you can use 
#' \code{parallel::detectCores()} to check how many 
#' cores are available on your computer
#' @param verbose if \code{TRUE} (default and recommended value), information
#' on the ongoing computations is printed in the console
#' 
#' @return A list containing the following objects:
#' \itemize{
#' \item \code{call}: the function call;
#' \item \code{concordance}: a data frame with the naive and
#' optimism-corrected estimates of the concordance (C) index;
#' \item \code{tdAUC}: a data frame with the naive and
#' optimism-corrected estimates of the time-dependent AUC
#' at the desired time points.
#' }
#' 
#' @import foreach doParallel glmnet survival survivalROC survcomp
#' @export
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
#' @seealso for the PRC-LMM model: \code{\link{fit_lmms}} (step 1),
#' \code{\link{summarize_lmms}} (step 2) and \code{\link{fit_prclmm}} (step 3);
#' for the PRC-MLPMM model: \code{\link{fit_mlpmms}} (step 1),
#' \code{\link{summarize_mlpmms}} (step 2) and \code{\link{fit_prcmlpmm}} (step 3).
#' 
#' @examples
#' \donttest{
# load fitted model
#' data(fitted_prclmm)
#' 
#' more.cores = FALSE
#' # IMPORTANT: set more.cores = TRUE to speed computations up!
#' if (!more.cores) n.cores = 2
#' if (more.cores) {
#'    # identify number of available cores on your machine
#'    n.cores = parallel::detectCores()
#'    if (is.na(n.cores)) n.cores = 2
#' }
#'                    
#' # compute the performance measures
#' perf = performance_prc(fitted_prclmm$step2, fitted_prclmm$step3, 
#'           times = c(0.5, 1, 1.5, 2), n.cores = n.cores)
#' 
#' # concordance index:
#' perf$concordance
#' # time-dependent AUC:
#' perf$tdAUC
#' }

performance_prc = function(step2, step3, times = 1,
                              n.cores = 1, verbose = TRUE) {
  call = match.call()
  # load namespaces
  requireNamespace('survival')
  requireNamespace('survcomp')
  requireNamespace('survivalROC')
  requireNamespace('glmnet')
  requireNamespace('foreach')
  # fix for 'no visible binding for global variable...' note
  i = b = NULL
  
  ############################
  ##### CHECK THE INPUTS #####
  ############################
  if (!is.numeric(times)) stop('times should be numeric!')
  n.times = length(times)
  c.out = data.frame(n.boots = NA, naive = NA,
                    cb.correction = NA, cb.performance = NA)
  tdauc.out = data.frame(pred.time = times, naive = NA,
                         cb.correction = NA, cb.performance = NA)
  # checks on step 2 input
  temp = c('call', 'ranef.orig', 'n.boots')
  check1 = temp %in% ls(step2)
  mess1 = paste('step2 input should cointain:', do.call(paste, as.list(temp)) )
  if (sum(check1) != 3) stop(mess1)
  n.boots = step2$n.boots
  ranef.orig = step2$ranef.orig
  # checks on step 3 input
  temp = c('call', 'pcox.orig', 'surv.data', 'n.boots')
  check2 = temp %in% ls(step3)
  mess2 = paste('step3 input should cointain:', do.call(paste, as.list(temp)) )
  if (sum(check2) != 4) stop(mess2)
  baseline.covs = step3$call$baseline.covs
  pcox.orig = step3$pcox.orig
  surv.data = step3$surv.data
  n = length(unique(surv.data$id))
  # further checks
  if (step2$n.boots != step3$n.boots) {
    stop('step2$n.boots and step3$n.boots are different!')
  }
  if (n.boots == 0) {
    mess = paste('The cluster bootstrap optimism correction has not',
            'been performed (n.boots = 0). Therefore, only the apparent',
            'values of the performance values will be returned.')
    warning(mess, immediate. = TRUE)
  }
  if (n.boots > 0) {
    temp = c('boot.ids', 'ranef.boot.train', 'ranef.boot.valid')
    check3 = temp %in% ls(step2)
    mess3 = paste('step2 should cointain:', do.call(paste, as.list(temp)) )
    if (sum(check3) != 3) stop(mess3)
    temp = c('boot.ids', 'pcox.boot')
    check4 = temp %in% ls(step3)
    mess4 = paste('step3 should cointain:', do.call(paste, as.list(temp)) )
    if (sum(check4) != 2) stop(mess4)
    boot.ids = step2$boot.ids
    ranef.boot.train = step2$ranef.boot.train
    ranef.boot.valid = step2$ranef.boot.valid
    pcox.boot = step3$pcox.boot
  }
  if (n.cores < 1) {
    warning('Input n.cores < 1, so we set n.cores = 1', immediate. = TRUE)
    n.cores = 1
  }
  # check how many cores are actually available for this computation
  if (n.boots > 0) {
    max.cores = parallel::detectCores()
    if (!is.na(max.cores)) {
      .check_ncores(avail = max.cores, requested = n.cores, verbose = verbose)
    }
  }
  
  # set up environment for parallel computing
  cl = parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)
  .info_ncores(n.cores, verbose = verbose)
  
  #############################
  ##### NAIVE PERFORMANCE #####
  #############################
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
  
  # C index on the original dataset
  relrisk.orig = predict(pcox.orig, 
                         newx = X.orig, 
                         s = 'lambda.min',
                         type = 'response')  
  c.naive = concordance.index(x = relrisk.orig, 
                    surv.time = surv.data$time,
                    surv.event = surv.data$event, 
                    method = "noether")
  c.out$n.boots = n.boots
  check = !inherits(c.naive, 'try-error')
  c.out$naive = ifelse (check, round(c.naive$c.index, 4), NA)
  
  # time-dependent AUC
  pmle.orig = as.numeric(coef(pcox.orig, s = 'lambda.min'))
  linpred.orig = X.orig %*% pmle.orig
  tdauc.naive = foreach(i = 1:n.times, .combine = 'c',
       .packages = c('survivalROC')) %dopar% {
    auc = try(survivalROC(Stime = surv.data$time, 
                          status = surv.data$event, 
                          marker = linpred.orig, 
                          entry = rep(0, n), 
                          predict.time = times[i], 
                          cut.values = NULL,
                          method = "NNE", 
                          span = 0.25*n^(-0.20)))
    check = !inherits(auc, 'try-error') & !is.nan(auc$AUC)
    out = ifelse(check, round(auc$AUC, 4), NA)
  }
  tdauc.out$naive = tdauc.naive

  ###############################
  ##### OPTIMISM CORRECTION #####
  ###############################
  if (n.boots > 0) {
    if (verbose) cat('Computation of optimism correction started\n')
    
    booty = foreach(b = 1:n.boots, .combine = 'rbind',
       .packages = c('survival', 'survcomp', 'survivalROC',
                     'glmnet', 'foreach')) %dopar% {
      # prepare data for the training set
      surv.data.train = surv.data[boot.ids[[b]], ]
      if (is.null(baseline.covs)) {
        X.train = as.matrix(ranef.boot.train[[b]])
      }
      if (!is.null(baseline.covs)) {
        X0 = model.matrix(as.formula(baseline.covs), data = surv.data.train)
        X.train = as.matrix(cbind(X0, ranef.boot.train[[b]]))
        contains.int = '(Intercept)' %in% colnames(X.train)
        if (contains.int) {
          X.train = X.train[ , -1] 
        }
      }
      # prepare X.valid for the validation set (surv.data already available)
      if (is.null(baseline.covs)) {
        X.valid = as.matrix(ranef.boot.valid[[b]])
      }
      if (!is.null(baseline.covs)) {
        X0 = model.matrix(as.formula(baseline.covs), data = surv.data)
        X.valid = as.matrix(cbind(X0, ranef.boot.valid[[b]]))
        contains.int = '(Intercept)' %in% colnames(X.valid)
        if (contains.int) {
          X.valid = X.valid[ , -1] 
        }
      }
      
      # C index on boot.train
      relrisk.train = predict(pcox.boot[[b]], 
                              newx = X.train, 
                              s = 'lambda.min',
                              type = 'response')  
      c.train = concordance.index(x = relrisk.train, 
                                  surv.time = surv.data.train$time,
                                  surv.event = surv.data.train$event, 
                                  method = "noether")
      
      # C index on boot.valid
      relrisk.valid = predict(pcox.boot[[b]], 
                              newx = X.valid, 
                              s = 'lambda.min',
                              type = 'response')  
      c.valid = concordance.index(x = relrisk.valid, 
                                  surv.time = surv.data$time,
                                  surv.event = surv.data$event, 
                                  method = "noether")
      
      # tdAUC on boot.train
      pmle.train = as.numeric(coef(pcox.boot[[b]], s = 'lambda.min'))
      linpred.train = X.train %*% pmle.train
      tdauc.train = foreach(i = 1:n.times, .combine = 'c') %do% {
        auc = try(survivalROC(Stime = surv.data.train$time, 
                              status = surv.data.train$event, 
                              marker = linpred.train, 
                              entry = rep(0, n), 
                              predict.time = times[i], 
                              cut.values = NULL,
                              method = "NNE", 
                              span = 0.25*n^(-0.20)))
        check = !inherits(auc, 'try-error') & !is.nan(auc$AUC)
        out = ifelse(check, round(auc$AUC, 4), NA)
      }
      
      # tdAUC on boot.valid
      pmle.train = as.numeric(coef(pcox.boot[[b]], s = 'lambda.min'))
      # important: the pmle comes from the model fitted on the training set
      linpred.valid = X.valid %*% pmle.train
      tdauc.valid = foreach(i = 1:n.times, .combine = 'c') %do% {
        auc = try(survivalROC(Stime = surv.data$time, 
                              status = surv.data$event, 
                              marker = linpred.valid, 
                              entry = rep(0, n), 
                              predict.time = times[i], 
                              cut.values = NULL,
                              method = "NNE", 
                              span = 0.25*n^(-0.20)))
        check = !inherits(auc, 'try-error') & !is.nan(auc$AUC)
        out = ifelse(check, round(auc$AUC, 4), NA)
      }
      
      # define outputs of parallel computing
      check1 = !inherits(c.train, 'try-error')
      ct = ifelse (check1, round(c.train$c.index, 4), NA)
      check2 = !inherits(c.valid, 'try-error')
      cv = ifelse (check2, round(c.valid$c.index, 4), NA)
      out = data.frame(repl = NA, stat = NA, times = NA, train = NA, valid = NA)
      out[1, ] = c(b, NA, NA, ct, cv)
      out[2:(n.times + 1),] = cbind(b, NA, times, tdauc.train, tdauc.valid)
      out$stat = c('C', rep('tdAUC', n.times))
      out$optimism = out$valid - out$train
      return(out)
    }
    
    # compute the optimism correction for the C index
    c.vals = booty[booty$stat == 'C', ]
    c.opt = mean(c.vals$optimism, na.rm = TRUE)
    c.out$cb.correction = round(c.opt, 4)
    c.out$cb.performance = c.out$naive + c.out$cb.correction
    
    # compute the optimism correction for the tdAUC
    tdauc.vals = booty[booty$stat == 'tdAUC', ]
    tdauc.opt = foreach(i = 1:n.times, .combine = 'c') %do% {
      temp = tdauc.vals[tdauc.vals$times == times[i], ]
      out = mean(temp$optimism, na.rm = TRUE)
      return(out)
    }
    tdauc.out$cb.correction = round(tdauc.opt, 4)
    tdauc.out$cb.performance = tdauc.out$naive + tdauc.out$cb.correction
    # closing message
    if (verbose) {
      cat('Computation of the optimism correction: finished :)\n')
    }
  }
  
  # close the cluster
  parallel::stopCluster(cl)

  names(c.out) = c('n.boots', 'C.naive', 'cb.opt.corr', 'C.adjusted')
  names(tdauc.out) = c('pred.time', 'tdAUC.naive', 'cb.opt.corr', 'tdAUC.adjusted')
  
  out = list('call' = call, 'concordance' = c.out, 
             'tdAUC' = tdauc.out)
  return(out)
}
