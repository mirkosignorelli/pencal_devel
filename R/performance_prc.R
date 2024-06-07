#' Predictive performance of the PRC-LMM and PRC-MLPMM models
#'
#' This function computes the naive and optimism-corrected
#' measures of performance (C index, time-dependent AUC and time-dependent 
#' Brier score) for the PRC models proposed 
#' in Signorelli et al. (2021). The optimism
#' correction is computed based on a cluster bootstrap
#' optimism correction procedure (CBOCP)
#' 
#' @param step2 the output of either \code{\link{summarize_lmms}} 
#' or \code{\link{summarize_mlpmms}} (step 2 of the estimation of
#' PRC)
#' @param step3 the output of \code{\link{fit_prclmm}} or
#' \code{\link{fit_prcmlpmm}} (step 3 of PRC)
#' @param metric the desired performance measure(s). Options include: 'tdauc',
#' 'c' and 'brier'
#' @param times numeric vector with the time points at which
#' to estimate the time-dependent AUC and time-dependent Brier score
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
#' at the desired time points;
#' \item \code{Brier}: a data frame with the naive and
#' optimism-corrected estimates of the time-dependent Brier score
#' at the desired time points;
#' }
#' 
#' @import foreach doParallel glmnet survival
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
#' # compute the time-dependent AUC
#' perf = performance_prc(fitted_prclmm$step2, fitted_prclmm$step3,
#'              metric = 'tdauc', times = c(3, 3.5, 4), n.cores = n.cores)
#'  # use metric = 'brier' for the Brier score and metric = 'c' for the
#'  # concordance index
#' 
#' # time-dependent AUC estimates:
#' ls(perf)
#' perf$tdAUC
#' }

performance_prc = function(step2, step3, metric = c('tdauc', 'c', 'brier'), 
                           times = c(2, 3), n.cores = 1, verbose = TRUE) {
  call = match.call()
  # load namespaces
  requireNamespace('survival')
  requireNamespace('glmnet')
  requireNamespace('foreach')
  # fix for 'no visible binding for global variable...' note
  i = b = model = NULL
  
  ############################
  ##### CHECK THE INPUTS #####
  ############################
  if (!is.numeric(times)) stop('times should be numeric!')
  n.times = length(times)
  compute.c = 'c' %in% metric
  compute.tdauc = 'tdauc' %in% metric
  compute.brier = 'brier' %in% metric
  c.out = data.frame(naive = NA, cb.correction = NA, cb.performance = NA)
  tdauc.out = data.frame(pred.time = times, naive = NA,
                         cb.correction = NA, cb.performance = NA)
  brier.out = data.frame(pred.time = times, naive = NA,
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
  if (max(times) >= max(surv.data$time)) {
    stop(paste('Some of the prediction times are larger than the larger event time in the dataset.',
                'Edit the times argument as appropriate'))
  }
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
  if (compute.c) {
    relrisk.orig = predict(pcox.orig, newx = X.orig, 
                           s = 'lambda.min', type = 'response') # exp(lin.pred)
    c.naive = survcomp::concordance.index(x = relrisk.orig, 
                                          surv.time = surv.data$time,
                                          surv.event = surv.data$event, 
                                          method = "noether")
    check = !inherits(c.naive, 'try-error')
    c.out$naive = ifelse (check, round(c.naive$c.index, 4), NA)
  }
  
  # time-dependent AUC
  pmle.orig = as.numeric(coef(pcox.orig, s = 'lambda.min'))
  linpred.orig = X.orig %*% pmle.orig
  if (compute.tdauc) {
    tdauc.naive = foreach(i = 1:n.times, .combine = 'c',
                          .packages = c('survivalROC')) %dopar% {
                            auc = try(survivalROC::survivalROC(Stime = surv.data$time, 
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
  }
  
  # time-dependent Brier score
  if (compute.brier) {
    # convert glmnet Cox model to equivalent with survival package
    df.orig = data.frame(time = surv.data$time,
                         event = surv.data$event,
                         linpred = linpred.orig,
                         id = surv.data$id)
    cox.survival = coxph(Surv(time = time, 
                              event = event) ~ linpred, 
                         data = df.orig, init = 1, 
                         control = coxph.control(iter.max = 0))
    # compute survival probabilities
    temp.sfit = survfit(cox.survival, newdata = df.orig,
                        se.fit = F, conf.int = F)
    spred = as.data.frame(t(summary(temp.sfit, times = times)$surv))
    # convert survival to failure probabilities and store them in a list
    fail_prob = as.list(1 - spred)
    # add names to the list
    names(fail_prob) = times
    # compute Brier Score and tdAUC
    perf = riskRegression::Score(fail_prob, times = times, metrics = 'brier',
                                 formula = Surv(time, event) ~ 1, data = surv.data,
                                 exact = FALSE, conf.int = FALSE, cens.model = "cox",
                                 splitMethod = "none", B = 0, verbose = FALSE)
    brier.naive = subset(perf$Brier$score, model == times)
    brier.out$naive = round(brier.naive$Brier, 4)
  }
  
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
      
      # C index
      if (compute.c) {
        # C index on boot.train
        relrisk.train = predict(pcox.boot[[b]], 
                                newx = X.train, 
                                s = 'lambda.min',
                                type = 'response')  
        c.train = survcomp::concordance.index(x = relrisk.train, 
                                              surv.time = surv.data.train$time,
                                              surv.event = surv.data.train$event, 
                                              method = "noether")
        # C index on boot.valid
        relrisk.valid = predict(pcox.boot[[b]], 
                                newx = X.valid, 
                                s = 'lambda.min',
                                type = 'response')  
        c.valid = survcomp::concordance.index(x = relrisk.valid, 
                                              surv.time = surv.data$time,
                                              surv.event = surv.data$event, 
                                              method = "noether")
      }

      #tdAUC
      pmle.train = as.numeric(coef(pcox.boot[[b]], s = 'lambda.min'))
      linpred.train = X.train %*% pmle.train
      linpred.valid = X.valid %*% pmle.train # important: the pmle comes from the model fitted on the training set

      if (compute.tdauc) {
        # tdAUC on boot.train
        tdauc.train = foreach(i = 1:n.times, .combine = 'c') %do% {
          auc = try(survivalROC::survivalROC(Stime = surv.data.train$time, 
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
        tdauc.valid = foreach(i = 1:n.times, .combine = 'c') %do% {
          auc = try(survivalROC::survivalROC(Stime = surv.data$time, 
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
      }
      
      # Brier score
      if (compute.brier) {
        # convert glmnet Cox model to equivalent with survival package
        df.train = data.frame(time = surv.data.train$time,
                             event = surv.data.train$event,
                             linpred = linpred.train,
                             id = surv.data.train$id)
        cox.train = coxph(Surv(time = time, 
                                  event = event) ~ linpred, 
                             data = df.train, init = 1, 
                             control = coxph.control(iter.max = 0))
        ### Brier on boot.train
        # compute survival probabilities
        temp.sfit = survfit(cox.train, newdata = df.train,
                            se.fit = F, conf.int = F)
        spred = as.data.frame(t(summary(temp.sfit, times = times)$surv))
        # convert survival to failure probabilities and store them in a list
        fail_prob.train = as.list(1 - spred)
        # add names to the list
        names(fail_prob.train) = times
        # compute Brier Score
        perf.train = riskRegression::Score(fail_prob.train, times = times, metrics = 'brier',
                                     formula = Surv(time, event) ~ 1, data = df.train,
                                     exact = FALSE, conf.int = FALSE, cens.model = "cox",
                                     splitMethod = "none", B = 0, verbose = FALSE)
        perf.train = subset(perf.train$Brier$score, model == times)
        brier.train = round(perf.train$Brier, 4)
        ### Brier on boot.valid
        # compute survival probabilities
        temp.sfit = survfit(cox.train, newdata = df.orig,
                            se.fit = F, conf.int = F)
        spred = as.data.frame(t(summary(temp.sfit, times = times)$surv))
        # convert survival to failure probabilities and store them in a list
        fail_prob.valid = as.list(1 - spred)
        # add names to the list
        names(fail_prob.valid) = times
        # compute Brier Score
        perf.valid = riskRegression::Score(fail_prob.valid, times = times, metrics = 'brier',
                                           formula = Surv(time, event) ~ 1, data = df.orig,
                                           exact = FALSE, conf.int = FALSE, cens.model = "cox",
                                           splitMethod = "none", B = 0, verbose = FALSE)
        perf.valid = subset(perf.valid$Brier$score, model == times)
        brier.valid = round(perf.valid$Brier, 4)
      }
      
      # define outputs of parallel computing
      out = data.frame(stat = NA, repl = NA, times = NA, train = NA, valid = NA)
      pos = 1
      if (compute.c) {
        check1 = !inherits(c.train, 'try-error')
        ct = ifelse (check1, round(c.train$c.index, 4), NA)
        check2 = !inherits(c.valid, 'try-error')
        cv = ifelse (check2, round(c.valid$c.index, 4), NA)
        out[pos, ] = c('C', b, NA, ct, cv)
        pos = pos + 1
      }
      if (compute.tdauc) {
        out[pos:(pos + n.times -1),] = cbind('tdAUC', b, times, tdauc.train, tdauc.valid)
        pos = pos + n.times
      }
      if (compute.brier) {
        out[pos:(pos + n.times -1),] = cbind('Brier', b, times, brier.train, brier.valid)
      }
      out[ , -1] = apply(out[ , -1], 2, as.numeric)
      out$optimism = out$valid - out$train
      return(out)
    }
    
    # compute the optimism correction for the C index
    if (compute.c) {
      c.vals = booty[booty$stat == 'C', ]
      c.opt = mean(c.vals$optimism, na.rm = TRUE)
      c.out$cb.correction = round(c.opt, 4)
      c.out$cb.performance = c.out$naive + c.out$cb.correction
    }
    
    # compute the optimism correction for the tdAUC
    if (compute.tdauc) {
      tdauc.vals = booty[booty$stat == 'tdAUC', ]
      tdauc.opt = foreach(i = 1:n.times, .combine = 'c') %do% {
        temp = tdauc.vals[tdauc.vals$times == times[i], ]
        out = mean(temp$optimism, na.rm = TRUE)
        return(out)
      }
      tdauc.out$cb.correction = round(tdauc.opt, 4)
      tdauc.out$cb.performance = tdauc.out$naive + tdauc.out$cb.correction
    }
    
    # compute the optimism correction for the Brier score
    if (compute.brier) {
      brier.vals = booty[booty$stat == 'Brier', ]
      brier.opt = foreach(i = 1:n.times, .combine = 'c') %do% {
        temp = brier.vals[brier.vals$times == times[i], ]
        out = mean(temp$optimism, na.rm = TRUE)
        return(out)
      }
      brier.out$cb.correction = round(brier.opt, 4)
      brier.out$cb.performance = brier.out$naive + brier.out$cb.correction
    }
    
    # closing message
    if (verbose) {
      cat('Computation of the optimism correction: finished :)\n')
    }
  }
  
  # close the cluster
  parallel::stopCluster(cl)
  
  # create outputs
  out = list('call' = call)
  if (compute.c) {
    names(c.out) = c('C.naive', 'optimism.correction', 'C.adjusted')
    out$concordance = c.out
  }
  if (compute.tdauc) {
    names(tdauc.out) = c('pred.time', 'tdAUC.naive', 'optimism.correction', 'tdAUC.adjusted')
    out$tdAUC = tdauc.out
  }
  if (compute.brier) {
    names(brier.out) = c('pred.time', 'Brier.naive', 'optimism.correction', 'Brier.adjusted')
    out$Brier = brier.out
  }
  return(out)
}
