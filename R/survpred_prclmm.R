#' Compute the predicted survival probabilities obtained
#' from the PRC models
#'
#' This function computes the predicted survival probabilities 
#' for the for the PRC-LMM model proposed 
#' in Signorelli et al. (2021)
#' 
#' @param step1 the output of \code{\link{fit_lmms}} (step 1
#' of the estimation of PRC-LMM)
#' @param step2 the output of \code{\link{summarize_lmms}} (step 2
#' of the estimation of PRC-LMM)
#' @param step3 the output of \code{\link{fit_prclmm}} (step 3
#' of the estimation of PRC-LMM)
#' @param times numeric vector with the time points at which
#' to estimate the time-dependent AUC
#' @param new.longdata longitudinal data if you want to compute 
#' predictions for new subjects on which the model was not trained.
#' It should comprise an identifier variable called `id`.
#' Default is \code{new.longdata = NULL}
#' @param new.basecovs a dataframe with baseline covariates for the
#' new subjects for which predictions are to be computed. 
#' It should comprise an identifier variable called `id`.
#'  Only needed if baseline covariates were included in step 3 and 
#' \code{new.longdata} is specified. Default is \code{new.basecovs = NULL}
#' @param keep.ranef should a data frame with the predicted random 
#' effects be included in the output? Default is \code{FALSE}
#' 
#' @return A list containing the function call (\code{call}),
#' a data frame with the predicted survival probabilities
#' computed at the supplied time points (\code{predicted_survival}),
#' and if \code{keep.ranef = TRUE} also the predicted random effects
#' \code{predicted_ranefs}.
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
#' @seealso \code{\link{fit_lmms}} (step 1),
#' \code{\link{summarize_lmms}} (step 2) and 
#' \code{\link{fit_prclmm}} (step 3)
#' 
#' @examples
#' # generate example data
#' set.seed(1234)
#' p = 4 # number of longitudinal predictors
#' simdata = simulate_prclmm_data(n = 100, p = p, p.relev = 2, 
#'              seed = 123, t.values = c(0, 0.2, 0.5, 1, 1.5, 2))
#'              
#' # step 1 of PRC-LMM: estimate the LMMs
#' y.names = paste('marker', 1:p, sep = '')
#' step1 = fit_lmms(y.names = y.names, 
#'                  fixefs = ~ age, ranefs = ~ age | id, 
#'                  long.data = simdata$long.data, 
#'                  surv.data = simdata$surv.data,
#'                  t.from.base = t.from.base,
#'                  n.boots = 0)
#'                  
#' # step 2 of PRC-LMM: compute the summaries 
#' # of the longitudinal outcomes
#' step2 = summarize_lmms(object = step1)
#' 
#' # step 3 of PRC-LMM: fit the penalized Cox models
#' step3 = fit_prclmm(object = step2, surv.data = simdata$surv.data,
#'                    baseline.covs = ~ baseline.age,
#'                    penalty = 'ridge')
#'                    
#' # predict survival probabilities at times 1, 2, 3
#' surv.probs = survpred_prclmm(step1, step2, step3, times = 1:3)
#' head(surv.probs$predicted_survival)
#' 
#' # predict survival probabilities for new subjects:
#' temp = simulate_prclmm_data(n = 10, p = p, p.relev = 2, 
#'       seed = 321, t.values = c(0, 0.2, 0.5, 1, 1.5, 2))
#' new.longdata = temp$long.data
#' new.basecovs = temp$surv.data[ , 1:2]
#' surv.probs.new = survpred_prclmm(step1, step2, step3, 
#'                      times = 1:3,
#'                      new.longdata = new.longdata,
#'                      new.basecovs = new.basecovs)
#' head(surv.probs.new$predicted_survival)

survpred_prclmm = function(step1, step2, step3, 
                           times = 1, new.longdata = NULL,
                           new.basecovs = NULL,
                           keep.ranef = FALSE) {
  call = match.call()
  # load namespaces
  requireNamespace('survival')
  requireNamespace('glmnet')
  requireNamespace('foreach')
  # fix for 'no visible binding for global variable...' note
  j = i = NULL
  
  ############################
  ##### CHECK THE INPUTS #####
  ############################
  if (!is.numeric(times)) stop('times should be numeric!')
  n.times = length(times)
  if (!is.logical(keep.ranef)) stop('keep.ranef should be T or F')
  
  # checks on step 1 input
  temp = c('call.info', 'lmm.fits.orig', 'df.sanitized', 'n.boots')
  check1 = temp %in% ls(step1)
  mess = paste('step1 input should cointain:',
               paste(temp, collapse = ', '))
  if (!all(check1, TRUE)) stop(mess)
  if (!is.null(new.longdata)) {
    new.longdata = new.longdata[order(new.longdata$id), ]
    vars = c(step1$call.info$y.names, 
             all.vars(step1$call.info$fixes),
             all.vars(step1$call.info$ranefs))
    check = all(vars %in% names(step1$df.sanitized), T)
    if (!check) {
      mess = 'new.longdata does not contain all variables employed in step 1'
      stop(mess)
    }
    if (!is.null(new.basecovs)) {
      new.basecovs = new.basecovs[order(new.basecovs$id), ]
    }
  }
  
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
  ### PRED RANEF WITH newdata ###
  ###############################
  if (!is.null(new.longdata)) {
    requireNamespace('nlme')
    # create matrix with predicted random effects computed on new dataset
    lmms = step1$lmm.fits.orig
    y.names = step1$call.info$y.names
    fixefs = step1$call.info$fixefs
    p = length(y.names)
    ids = unique(new.longdata$id)
    n = length(ids)
    
    ranef.new = foreach(j = 1:p, .combine = 'cbind') %do% {
      # check if NAs on response
      new.df = new.longdata
      y = new.df[ , y.names[j]]
      nas = which(is.na(y))
      if (length(nas) > 0) {
        new.df = new.df[-nas, ]
        y = new.df[ , y.names[j]]
      }
      # check that we didn't lose subjects
      check = (length(unique(new.df$id)) == n)
      if (!check) {
        lost_ids = setdiff(new.longdata$id, new.df$id)
        mess = paste('Variable ', y.names[j], ': all values are NA for at least 1 subject.', 
                     'Predictions obtained by setting predicted random effects for such subjects = 0 (population average)',
                     sep = '')
        warning(mess)
      }
      # retrieve the right pieces from lmm
      D.hat = getVarCov(lmms[[j]], type = 'random.effects')
      beta.hat = fixef(lmms[[j]])
      sigma2.hat = summary(lmms[[j]])$sigma^2 # error variance
      # create X and Z
      X = model.matrix(as.formula(fixefs), data = new.df)
      formYz = formula(lmms[[j]]$modelStruct$reStruct[[1]])
      mfZ = model.frame(terms(formYz), data = new.df)
      Z = model.matrix(formYz, mfZ)
      #####
      current = foreach(i = 1:n, .combine = 'rbind') %do% {
        rows = which(new.df$id == ids[i])
        if (length(rows) > 0) {
          Xi = X[rows, , drop = FALSE] 
          # drop = F prevents conversion to vector when rows has length 1
          yi = y[rows]
          Zi = Z[rows, , drop = FALSE]
          I.matr = diag(1, length(rows), length(rows))
          Vi = Zi %*% D.hat %*% t(Zi) + sigma2.hat * I.matr
          temp = yi - Xi %*% beta.hat
          t(D.hat %*% t(Zi) %*% solve(Vi) %*% temp)
        }
        else {
          matrix(0, 1, ncol(Z))
        }
      }
      # return
      current
    }
    rownames(ranef.new) = ids
  }
  
  ###############################
  ##### COMPUTE PREDICTIONS #####
  ###############################
  # reconstruct pieces  
  if (is.null(baseline.covs)) {
    X.orig = as.matrix(ranef.orig)
    rownames(X.orig) = rownames(ranef.orig) # fixed in v 1.0.3
    if (!is.null(new.longdata)) {
      X.new = as.matrix(ranef.new)
      rownames(X.new) = ids # fixed in v 1.0.3
    }
  }
  if (!is.null(baseline.covs)) {
    X0 = model.matrix(as.formula(baseline.covs), data = surv.data)
    X.orig = as.matrix(cbind(X0, ranef.orig))
    rownames(X.orig) = rownames(ranef.orig) # fixed in v 1.0.3
    contains.int = '(Intercept)' %in% colnames(X.orig)
    if (contains.int) {
      X.orig = X.orig[ , -1, drop = FALSE] 
    }
    if (!is.null(new.longdata)) {
      if (is.null(new.basecovs)) {
        stop('new.basecovs was not supplied')
      }
      X0 = model.matrix(as.formula(baseline.covs), 
                        data = new.basecovs)
      X.new = as.matrix(cbind(X0, ranef.new))
      rownames(X.new) = ids # fixed in v 1.0.3
      contains.int = '(Intercept)' %in% colnames(X.new)
      if (contains.int) {
        X.new = X.new[ , -1, drop = FALSE]
      }
    }
  }
  
  beta.hat = coef(pcox.orig, s = 'lambda.min')
  temp = rownames(beta.hat)
  beta.hat = as.numeric(beta.hat)
  names(beta.hat) = temp
  # check if dimensionality correct
  if (ncol(X.orig) != length(beta.hat)) stop("Dimensions of beta.hat and X.train don't correspond")

  # convert glmnet Cox model to equivalent with survival package
  df.orig = data.frame(time = surv.data$time,
                       event = surv.data$event,
                       linpred = X.orig %*% beta.hat,
                       id = surv.data$id)
  cox.survival = coxph(Surv(time = time, 
                            event = event) ~ linpred, 
                       data = df.orig, init = 1, 
     control = coxph.control(iter.max = 0))
  
  # compute survival probabilities
  if (is.null(new.longdata)) X.new = X.orig
  df.new = data.frame(linpred = X.new %*% beta.hat)
  temp.sfit = survfit(cox.survival,
       newdata = df.new, se.fit = F, conf.int = F)
  shat.orig = t(summary(temp.sfit, 
                       times = times)$surv)
  # output
  preds = data.frame(rownames(X.new), shat.orig)
  names(preds) = c('id', paste('S(', times, ')', sep = '') )
  out = list('call' = call, 
             'predicted_survival' = preds)
  if (keep.ranef) {
    pred.ranef = ranef.orig
    if (!is.null(new.longdata)) {
      pred.ranef = as.data.frame(ranef.new)
      names(pred.ranef) = colnames(ranef.orig)
    }
    out = list('call' = call,
               'predicted_survival' = preds,
               'predicted_ranefs' = pred.ranef)
  }
  return(out)
}
