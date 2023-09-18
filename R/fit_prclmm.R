#' Step 3 of PRC-LMM (estimation of the penalized Cox model(s))
#'
#' This function performs the third step for the estimation
#' of the PRC-LMM model proposed in Signorelli et al. (2021)
#' 
#' @param object the output of step 2 of the PRC-LMM procedure, 
#' as produced by the \code{\link{summarize_lmms}} function
#' @param surv.data a data frame with the survival data and (if 
#' relevant) additional baseline covariates. \code{surv.data} should at least
#' contain a subject id (called \code{id}), the time to event outcome  
#' (\code{time}), and binary event variable (\code{event})
#' @param baseline.covs a formula specifying the variables 
#' (e.g., baseline age) in \code{surv.data} that should be included 
#' as baseline covariates in the penalized Cox model. Example:
#' \code{baseline.covs = '~ baseline.age'}. Default is \code{NULL}
#' @param penalty the type of penalty function used for regularization.
#' Default is \code{'ridge'}, other possible values are \code{'elasticnet'} 
#' and \code{'lasso'}
#' @param standardize logical argument: should the predictors (both baseline covariates
#' and predicted random effects) be standardized when included  as covariates
#' in the penalized Cox model? Default is \code{TRUE}
#' @param pfac.base.covs a single value, or a vector of values, indicating
#' whether the baseline covariates (if any) should be penalized (1) or not (0).
#' Default is \code{pfac.base.covs = 0} (no penalization of all baseline covariates)
#' @param n.alpha.elnet number of alpha values for the two-dimensional 
#' grid of tuning parameteres in elasticnet.
#' Only relevant if \code{penalty = 'elasticnet'}. Default is 11,
#' so that the resulting alpha grid is c(1, 0.9, 0.8, ..., 0.1, 0)
#' @param n.folds.elnet number of folds to be used for the selection
#' of the tuning parameter in elasticnet. Only relevant if 
#' \code{penalty = 'elasticnet'}. Default is 5
#' @param n.cores number of cores to use to parallelize part of
#' the computations. If \code{ncores = 1} (default), no parallelization is done. 
#' Pro tip: you can use \code{parallel::detectCores()} to check 
#' how many cores are available on your computer
#' @param verbose if \code{TRUE} (default and recommended value), information
#' on the ongoing computations is printed in the console
#' 
#' @return A list containing the following objects:
#' \itemize{
#' \item \code{call}: the function call
#' \item \code{pcox.orig}: the penalized Cox model fitted on the
#' original dataset;
#' \item \code{surv.data}: the supplied survival data (ordered by
#' subject id)
#' \item \code{n.boots}: number of bootstrap samples;
#' \item \code{boot.ids}: a list with the ids of bootstrapped subjects 
#' (when \code{n.boots > 0});
#' \item \code{pcox.boot}: a list where each element is a fitted penalized
#' Cox model for a given bootstrap sample (when \code{n.boots > 0}).
#' }
#' 
#' @import foreach doParallel survival glmnet stats
#' @importFrom dplyr arrange
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
#' \code{\link{summarize_lmms}} (step 2),
#' \code{\link{performance_prc}}
#' 
#' @examples
#' # generate example data
#' set.seed(1234)
#' p = 4 # number of longitudinal predictors
#' simdata = simulate_prclmm_data(n = 100, p = p, p.relev = 2, 
#'              seed = 123, t.values = c(0, 0.2, 0.5, 1, 1.5, 2))
#'              
#' # specify options for cluster bootstrap optimism correction
#' # procedure and for parallel computing 
#' do.bootstrap = FALSE
#' # IMPORTANT: set do.bootstrap = TRUE to compute the optimism correction!
#' n.boots = ifelse(do.bootstrap, 100, 0)
#' more.cores = FALSE
#' # IMPORTANT: set more.cores = TRUE to parallelize and speed computations up!
#' if (!more.cores) n.cores = 1
#' if (more.cores) {
#'    # identify number of available cores on your machine
#'    n.cores = parallel::detectCores()
#'    if (is.na(n.cores)) n.cores = 8
#' }
#' 
#' # step 1 of PRC-LMM: estimate the LMMs
#' y.names = paste('marker', 1:p, sep = '')
#' step1 = fit_lmms(y.names = y.names, 
#'                  fixefs = ~ age, ranefs = ~ age | id, 
#'                  long.data = simdata$long.data, 
#'                  surv.data = simdata$surv.data,
#'                  t.from.base = t.from.base,
#'                  n.boots = n.boots, n.cores = n.cores)
#'                  
#' # step 2 of PRC-LMM: compute the summaries 
#' # of the longitudinal outcomes
#' step2 = summarize_lmms(object = step1, n.cores = n.cores)
#' 
#' # step 3 of PRC-LMM: fit the penalized Cox models
#' step3 = fit_prclmm(object = step2, surv.data = simdata$surv.data,
#'                    baseline.covs = ~ baseline.age,
#'                    penalty = 'ridge', n.cores = n.cores)
#' summary(step3)                    

fit_prclmm = function(object, surv.data, baseline.covs = NULL,
                      penalty = 'ridge', standardize = TRUE,
                      pfac.base.covs = 0,
                      n.alpha.elnet = 11, n.folds.elnet = 5,
                      n.cores = 1, verbose = TRUE) {
  call = match.call()
  requireNamespace('foreach')
  requireNamespace('glmnet')
  requireNamespace('survival')
  requireNamespace('doParallel')
  # fix for 'no visible binding for global variable...' note
  id = i = b = NULL
  # identify inputs and perform checks
  ranef.orig = object$ranef.orig
  n.boots = object$n.boots
  do.bootstrap = ifelse(n.boots > 0, TRUE, FALSE)
  if (do.bootstrap) {
    boot.ids = object$boot.ids
    ranef.boot.train = object$ranef.boot.train
  }
  # preliminary checks
  check1 = is.data.frame(surv.data)
  if (!check1) stop('surv.data should be a data.frame')
  check2 = c('id', 'time', 'event') %in% names(surv.data)
  if (sum(check2) != 3) stop("surv.data should contain at least: 'id', 'time', 'event' variables")
  # order survdata by id (in case it wasn't sorted) because
  # in step 1 everything was also sorted by id (relevant for bootstrap!)
  surv.data = dplyr::arrange(surv.data, id)
  check3 = penalty %in% c('ridge', 'elasticnet', 'lasso')
  if (!check3) stop('penalty should be ridge, elasticnet or lasso')
  if (penalty == 'elasticnet') {
    if (n.alpha.elnet < 3) {
      stop('Supply n.alpha.elnet >= 3 so that a two-dimensional (alpha, lambda) grid 
           to select the elasticnet tuning parameters can be created')  
    }
    if (n.alpha.elnet > 21) {
      warning('n.alpha.elnet > 21 may significantly slow computations down.
            Consider setting n.alpha.elnet to a smaller value', 
              immediate. = TRUE)  
    }
    if (n.folds.elnet < 5) {
      stop('Doing cross-validation with less than 5 folds is not recommended') 
    }
    if (n.folds.elnet > 10) {
      warning('Whoa! You are using a lot of folds to select the elasticnet tuning
              parameters! Are you sure you want to do that? Please proceed 
              with care', immediate. = TRUE) 
    }
  }
  check8 = is.logical(standardize)
  if (!check8) stop('standardize must be a logical argument (T / F)')
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
  
  # prepare environment depending on which penalty is chosen
  double.cv = ifelse(penalty == 'elasticnet', TRUE, FALSE)
  if (!double.cv) { # ridge or lasso
    alpha = ifelse(penalty == 'lasso', 1, 0)
  }
  if (double.cv) { # elasticnet
    alpha = seq(1, 0, length.out = n.alpha.elnet)
    n.reps = floor(nrow(ranef.orig) / n.folds.elnet)
    fold.ids = rep(1:n.folds.elnet, n.reps)
    n.remaining = nrow(ranef.orig) - n.folds.elnet*n.reps
    if (n.remaining > 0) {
      fold.ids = c(fold.ids, 1:n.remaining)
    }
  }
  
  ###########################################
  ###### fit step 3 on original dataset #####
  ###########################################
  if (verbose) cat('Estimating penalized Cox model on the original dataset...\n')
  surv.orig = Surv(time = surv.data$time, event = surv.data$event)
  if (is.null(baseline.covs)) {
    X.orig = as.matrix(ranef.orig)
    pen.fac = rep(1, ncol(X.orig))
  }
  if (!is.null(baseline.covs)) {
    X0 = model.matrix(as.formula(baseline.covs), data = surv.data)
    contains.int = '(Intercept)' %in% colnames(X0)
    if (contains.int) {
      X0 = X0[ , -1, drop = FALSE] 
    }
    X.orig = as.matrix(cbind(X0, ranef.orig))
    # fix penalty factor
    if (length(pfac.base.covs) > ncol(X0)) {
      stop('pfac.base.covs contains too many elements!')
    }
    if (length(pfac.base.covs) < ncol(X0) & length(pfac.base.covs) > 1) {
      stop('wrong number of elements in pfac.base.covs')
    }
    if (length(pfac.base.covs) == 1) {
      pfac.base.covs = rep(pfac.base.covs, ncol(X0)-1)
    }
    pen.fac = c(pfac.base.covs, rep(1, ncol(X.orig) - length(pfac.base.covs)))
  }
  
  # if penalty is ridge or lasso:
  if (!double.cv) {
    # alpha already defined above :)
    pcox.orig = cv.glmnet(x = X.orig, y = surv.orig, family="cox", 
                     standardize = standardize, 
                     penalty.factor = pen.fac,
                     alpha = alpha, maxit = 1e7)
  }
  # if penalty is elasticnet:
  if (double.cv) {
    fits = vector('list', n.alpha.elnet)
    tuning.matr = matrix(NA, n.alpha.elnet, 3)
    colnames(tuning.matr) = c('alpha', 'lambda.min', 'deviance')
    # run cross-validation
    foreach (i = 1:n.alpha.elnet) %do% {
      alpha.el = alpha[i]
      temp = cv.glmnet(x = X.orig, y = surv.orig, family="cox", 
                            standardize = standardize,
                            penalty.factor = pen.fac,
                            foldid = fold.ids,
                            alpha = alpha.el, maxit = 1e7)
      tuning.matr[i, ] = c(alpha.el, temp$lambda.min,
                           min(temp$cvm))
      fits[[i]] = temp
    }
    # optimal combination of (alpha, lambda)
    id.best = which.min(tuning.matr[ ,3])
    pcox.orig = fits[[id.best]]
  }
  if (verbose) cat('...done\n')
  
  #######################
  ### start bootstrap ###
  #######################
  if (do.bootstrap) {
    if (verbose) cat('Bootstrap procedure started\n')
    # set up environment for parallel computing
    cl = parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    .info_ncores(n.cores, verbose = verbose)
    
    booty = foreach(b = 1:n.boots,
    .packages = c('survival', 'pencal', 'glmnet')) %dopar% {
      # prepare data
      boot.surv.data = surv.data[boot.ids[[b]], ]
      surv.boot = Surv(time = boot.surv.data$time, 
                       event = boot.surv.data$event)
      if (is.null(baseline.covs)) {
        X.boot = as.matrix(ranef.boot.train[[b]])
        pen.fac = rep(1, ncol(X.boot))
      }
      if (!is.null(baseline.covs)) {
        X0 = model.matrix(as.formula(baseline.covs), data = boot.surv.data)
        X.boot = as.matrix(cbind(X0, ranef.boot.train[[b]]))
        contains.int = '(Intercept)' %in% colnames(X.boot)
        if (contains.int) {
          X.boot = X.boot[ , -1] 
        }
        pen.fac = c(pfac.base.covs, rep(1, ncol(X.boot) - length(pfac.base.covs)))
      }

      # if penalty is ridge or lasso:
      if (!double.cv) {
        # alpha already defined above :)
        pcox.boot = cv.glmnet(x = X.boot, y = surv.boot, family="cox", 
                              standardize = standardize,
                              penalty.factor = pen.fac,
                              alpha = alpha, maxit = 1e7)
      }
      # if penalty is elasticnet:
      if (double.cv) {
        fits = vector('list', n.alpha.elnet)
        tuning.matr = matrix(NA, n.alpha.elnet, 3)
        colnames(tuning.matr) = c('alpha', 'lambda.min', 'deviance')
        # run cross-validation
        foreach (i = 1:n.alpha.elnet) %do% {
          alpha.el = alpha[i]
          temp = cv.glmnet(x = X.boot, y = surv.boot, family="cox", 
                           standardize = standardize,
                           penalty.factor = pen.fac,
                           foldid = fold.ids,
                           alpha = alpha.el, maxit = 1e7)
          tuning.matr[i, ] = c(alpha.el, temp$lambda.min,
                               min(temp$cvm))
          fits[[i]] = temp
        }
        # optimal combination of (alpha, lambda)
        id.best = which.min(tuning.matr[ ,3])
        pcox.boot = fits[[id.best]]
      }
      # return statement for foreach:
      out = list('pcox.boot' = pcox.boot, 
                 'surv.boot' = surv.boot, 
                 'X.boot' = X.boot)
      out
    }
    # re-format exports
    pcox.boot = foreach (b = 1:n.boots) %dopar% {
      booty[[b]]$pcox.boot
    }
    surv.boot = foreach (b = 1:n.boots) %dopar% {
      booty[[b]]$surv.boot
    }
    X.boot = foreach (b = 1:n.boots) %dopar% {
      booty[[b]]$X.boot
    }
    # close the cluster
    parallel::stopCluster(cl)  
    if (verbose) {
      cat('Bootstrap procedure finished\n')
      cat('Computation of step 3: finished :)\n')
    }
  }
  
  # export results
  out = list('call' = call, 'pcox.orig' = pcox.orig,
            'surv.data' = surv.data, 'n.boots' = n.boots)
  if (n.boots >= 1) {
    out[['boot.ids']] = boot.ids
    out[['pcox.boot']] = pcox.boot
  }
  class(out) = 'prclmm'
  return(out)
}
