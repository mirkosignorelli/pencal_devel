#' Predictive performance of the penalized Cox model
#' with time-independent covariates
#'
#' This function computes the naive and optimism-corrected
#' measures of performance (C index, time-dependent AUC and time-dependent 
#' Brier score) for a penalized Cox model with time-independent covariates. 
#' The optimism correction is computed based on a cluster bootstrap
#' optimism correction procedure (CBOCP, Signorelli et al., 2021)
#' 
#' @param fitted_pencox the output of \code{\link{pencox}}
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
#' at the desired time points.
#' }
#' 
#' @import foreach doParallel glmnet survival
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
#' @seealso \code{\link{pencox}}
#' 
#' @examples
#' # generate example data
#' set.seed(1234)
#' p = 4 # number of longitudinal predictors
#' simdata = simulate_prclmm_data(n = 100, p = p, p.relev = 2, 
#'              seed = 123, t.values = c(0, 0.2, 0.5, 1, 1.5, 2))
#' # create dataframe with baseline measurements only
#' baseline.visits = simdata$long.data[which(!duplicated(simdata$long.data$id)),]
#' df = cbind(simdata$surv.data, baseline.visits)
#' df = df[ , -c(5:7)]
#' 
#' do.bootstrap = FALSE
#' # IMPORTANT: set do.bootstrap = TRUE to compute the optimism correction!
#' n.boots = ifelse(do.bootstrap, 100, 0)
#' more.cores = FALSE
#' # IMPORTANT: set more.cores = TRUE to speed computations up!
#' if (!more.cores) n.cores = 2
#' if (more.cores) {
#'    # identify number of available cores on your machine
#'    n.cores = parallel::detectCores()
#'    if (is.na(n.cores)) n.cores = 2
#' }
#' 
#' form = as.formula(~ baseline.age + marker1 + marker2
#'                      + marker3 + marker4)
#' base.pcox = pencox(data = df, 
#'               formula = form, 
#'               n.boots = n.boots, n.cores = n.cores) 
#' ls(base.pcox)
#'                    
#' # compute the performance measures
#' perf = performance_pencox(fitted_pencox = base.pcox, 
#'           metric = 'tdauc', times = c(1, 1.5, 2), n.cores = n.cores)
#'  # use metric = 'brier' for the Brier score and metric = 'c' for the
#'  # concordance index
#' 
#' # time-dependent AUC estimates:
#' ls(perf)
#' perf$tdAUC

performance_pencox = function(fitted_pencox, metric = c('tdauc', 'c', 'brier'), 
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
  # checks on fitted_pencox
  temp = c('call', 'pcox.orig', 'surv.data', 'X.orig', 'n.boots')
  check2 = temp %in% ls(fitted_pencox)
  mess2 = paste('fitted_pencox input should cointain:', 
                do.call(paste, as.list(temp)) )
  if (sum(check2) != 5) stop(mess2)
  X.orig = fitted_pencox$X.orig
  pcox.orig = fitted_pencox$pcox.orig
  surv.data = fitted_pencox$surv.data
  n = length(unique(surv.data$id))
  n.boots = fitted_pencox$n.boots
  # further checks
  if (n.boots == 0) {
    mess = paste('The bootstrap optimism correction has not',
            'been performed (n.boots = 0). Therefore, only the apparent',
            'values of the performance values will be returned.')
    warning(mess, immediate. = TRUE)
  }
  if (n.boots > 0) {
    temp = c('boot.ids', 'pcox.boot')
    check3 = temp %in% ls(fitted_pencox)
    mess3 = paste('fitted_pencox should cointain:', do.call(paste, as.list(temp)) )
    if (sum(check3) != 2) stop(mess3)
    boot.ids = fitted_pencox$boot.ids
    pcox.boot = fitted_pencox$pcox.boot
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
  
  # C index on the original dataset
  if (compute.c) {
    relrisk.orig = predict(pcox.orig, 
                           newx = X.orig, 
                           s = 'lambda.min',
                           type = 'response')  
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
      X.train = X.orig[boot.ids[[b]], ]
      X.valid = X.orig
      
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
      
      pmle.train = as.numeric(coef(pcox.boot[[b]], s = 'lambda.min'))
      linpred.train = X.train %*% pmle.train
      linpred.valid = X.valid %*% pmle.train
      
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
