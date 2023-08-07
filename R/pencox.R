#' Estimation of a penalized Cox model with time-independent covariates
#'
#' This function estimates a penalized Cox model where only
#' time-independent covariates are included as predictors, and then
#' computes a bootstrap optimism correction procedure that 
#' is used to validate the predictive performance of the model
#' 
#' @param data a data frame with one row for each subject.It
#' should at least contain a subject id (called \code{id}), 
#' the time to event outcome (\code{time}), and the binary censoring
#' indicator (\code{event}), plus at least one covariate to
#' be included in the linear predictor
#' @param formula a formula specifying the variables 
#' in \code{data} to include as predictors in 
#' the penalized Cox model
#' @param penalty the type of penalty function used for regularization.
#' Default is \code{'ridge'}, other possible values are \code{'elasticnet'} 
#' and \code{'lasso'}
#' @param standardize logical argument: should the covariates
#' be standardized when included in the penalized Cox model? Default is \code{TRUE}
#' @param penalty.factor a single value, or a vector of values, indicating
#' whether the covariates (if any) should be penalized (1) or not (0).
#' Default is \code{penalty.factor = 1}
#' @param n.alpha.elnet number of alpha values for the two-dimensional 
#' grid of tuning parameteres in elasticnet.
#' Only relevant if \code{penalty = 'elasticnet'}. Default is 11,
#' so that the resulting alpha grid is c(1, 0.9, 0.8, ..., 0.1, 0)
#' @param n.folds.elnet number of folds to be used for the selection
#' of the tuning parameter in elasticnet. Only relevant if 
#' \code{penalty = 'elasticnet'}. Default is 5
#' @param n.boots number of bootstrap samples to be used 
#' in the bootstrap optimism correction procedure. If 0, no
#' bootstrapping is performed
#' @param n.cores number of cores to use to parallelize the computation
#' of the CBOCP. If \code{ncores = 1} (default), no parallelization is done. 
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
#' \item \code{surv.data}: a data frame with the survival data
#' \item \code{X.orig}: a data frame with the design matrix used
#' to estimate the Cox model
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
#' @seealso \code{\link{fit_prclmm}}, 
#' \code{\link{fit_prcmlpmm}}
#' 
#' @examples
#' # generate example data
#' set.seed(1234)
#' p = 4 # number of longitudinal predictors
#' simdata = simulate_prclmm_data(n = 100, p = p, p.relev = 2, 
#'              seed = 123, t.values = c(0, 0.2, 0.5, 1, 1.5, 2))
#' #create dataframe with baseline measurements only
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

pencox = function(data, formula,
                      penalty = 'ridge', standardize = TRUE,
                      penalty.factor = 1,
                      n.alpha.elnet = 11, n.folds.elnet = 5,
                      n.boots = 0, n.cores = 1, verbose = TRUE) {
  call = match.call()
  requireNamespace('foreach')
  requireNamespace('glmnet')
  requireNamespace('survival')
  requireNamespace('doParallel')
  # fix for 'no visible binding for global variable...' note
  id = i = b = NULL
  # preliminary checks
  check1 = is.data.frame(data)
  if (!check1) stop('data should be a data.frame')
  check2 = c('id', 'time', 'event') %in% names(data)
  if (sum(check2) != 3) stop("data should contain at least: 'id', 'time', 'event' variables")
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
  check9 = inherits(formula, 'formula')
  if (!check9) stop('formula should be of class: formula')
  
  # check how many cores are actually available for this computation
  do.bootstrap = ifelse(n.boots > 0, TRUE, FALSE)
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
    n.reps = floor(nrow(data) / n.folds.elnet)
    fold.ids = rep(1:n.folds.elnet, n.reps)
    n.remaining = nrow(data) - n.folds.elnet*n.reps
    if (n.remaining > 0) {
      fold.ids = c(fold.ids, 1:n.remaining)
    }
  }
  
  ###########################################
  ###### fit penCox on original dataset #####
  ###########################################
  if (verbose) cat('Estimated penalized Cox model on the original dataset...\n')
  surv.orig = Surv(time = data$time, event = data$event)
  X.orig = model.matrix(formula, data)[ , -1]
  # fix penalty factor
  if (length(penalty.factor) > ncol(X.orig)) {
    stop('penalty.factor contains too many elements!')
  }
  if (length(penalty.factor) < ncol(X.orig) & length(penalty.factor) > 1) {
    stop('wrong number of elements in penalty.factor')
  }
  if (length(penalty.factor) == 1) {
    penalty.factor = rep(penalty.factor, ncol(X.orig))
  }

  # if penalty is ridge or lasso:
  if (!double.cv) {
    # alpha already defined above :)
    pcox.orig = cv.glmnet(x = X.orig, y = surv.orig, family="cox", 
                     standardize = standardize, 
                     penalty.factor = penalty.factor,
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
                            penalty.factor = penalty.factor,
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
    
    # sample bootstrap row ids
    n = nrow(data)
    boot.ids = foreach(i = 1:n.boots) %do% {
      set.seed(i)
      sort(sample(1:n, n, TRUE))
    }
    
    booty = foreach(b = 1:n.boots,
    .packages = c('survival', 'pencal', 'glmnet')) %dopar% {
      # prepare data
      boot.surv.data = data[boot.ids[[b]], ]
      surv.boot = Surv(time = boot.surv.data$time, 
                       event = boot.surv.data$event)
      X.boot = model.matrix(formula, 
                            data = boot.surv.data)[ , -1] 

      # if penalty is ridge or lasso:
      if (!double.cv) {
        # alpha already defined above :)
        pcox.boot = cv.glmnet(x = X.boot, y = surv.boot, family="cox", 
                              standardize = standardize,
                              penalty.factor = penalty.factor,
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
                           penalty.factor = penalty.factor,
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
      cat('Computation of bootstrap procedure: finished :)\n')
    }
  }
  
  # export results
  surv.data = data[, c('id', 'time', 'event')]
  out = list('call' = call, 'pcox.orig' = pcox.orig,
            'surv.data' = surv.data, 'X.orig' = X.orig,
            'n.boots' = n.boots)
  if (n.boots >= 1) {
    out[['boot.ids']] = boot.ids
    out[['pcox.boot']] = pcox.boot
  }
  return(out)
}
