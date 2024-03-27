#' Step 2 of PRC-MLPMM (computation of the predicted random effects)
#'
#' This function performs the second step for the estimation
#' of the PRC-MLPMM model proposed in Signorelli et al. (2021)
#' 
#' @param object a list of objects as produced by \code{\link{fit_mlpmms}}
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
#' \item \code{ranef.orig}: a matrix with the predicted random effects
#' computed for the original data;
#' \item \code{n.boots}: number of bootstrap samples;
#' \item \code{boot.ids}: a list with the ids of bootstrapped subjects 
#' (when \code{n.boots > 0});
#' \item \code{ranef.boot.train}: a list where each element is a matrix that 
#' contains the predicted random effects for each bootstrap sample 
#' (when \code{n.boots > 0});
#' \item \code{ranef.boot.valid}: a list where each element is a matrix that 
#' contains the predicted random effects on the original data, based on the 
#' mlpmms fitted on the cluster bootstrap samples (when \code{n.boots > 0});
#' }
#' 
#' @import foreach doParallel stats
#' @export
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
#' 
#' @seealso \code{\link{fit_mlpmms}} (step 1), 
#' \code{\link{fit_prcmlpmm}} (step 3),
#' \code{\link{performance_prc}}
#' 
#' @examples
#' \donttest{
#' # generate example data
#' set.seed(123)
#' n.items = c(4,2,2,3,4,2)
#' simdata = simulate_prcmlpmm_data(n = 100, p = length(n.items),  
#'              p.relev = 3, n.items = n.items, 
#'              type = 'u+b', seed = 1)
#'  
#' # specify options for cluster bootstrap optimism correction
#' # procedure and for parallel computing 
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
#' # step 1 of PRC-MLPMM: estimate the MLPMMs
#' y.names = vector('list', length(n.items))
#' for (i in 1:length(n.items)) {
#'   y.names[[i]] = paste('marker', i, '_', 1:n.items[i], sep = '')
#' }
#' 
#' step1 = fit_mlpmms(y.names, fixefs = ~ contrast(age),  
#'                  ranef.time = age, randint.items = TRUE, 
#'                  long.data = simdata$long.data, 
#'                  surv.data = simdata$surv.data,
#'                  t.from.base = t.from.base,
#'                  n.boots = n.boots, n.cores = n.cores)
#'
#' # step 2 of PRC-MLPMM: compute the summaries 
#' step2 = summarize_mlpmms(object = step1, n.cores = n.cores)
#' summary(step2)
#' }

summarize_mlpmms = function(object, n.cores = 1, verbose = TRUE) {
  # fix for 'no visible binding for global variable...' note
  j = NULL
  
  call = match.call()
  # load namespaces
  requireNamespace('lcmm')
  requireNamespace('foreach')
  requireNamespace('doParallel')
  
  # checks on step 1 output + retrieve info from it
  if (!is.list(object)) stop('object should be a list, as outputted from fit_lmms')
  minimal.elements = c('call.info', 'mlpmm.fits.orig', 'df.sanitized', 'n.boots')
  check1 = minimal.elements %in% ls(object)
  mess = paste('At least one of the following elements is missing in object:',
               paste(minimal.elements, collapse = ', '))
  if (!all(check1, TRUE)) stop(mess)
  df.sanitized = object$df.sanitized
  y.names = object$call.info$y.names
  fixefs = object$call.info$fixefs
  ranefs = object$call.info$ranefs
  randint.items = object$call.info$randint.items
  p = length(y.names)
  n.boots = object$n.boots
  if (n.boots < 0) {
    warning('Input n.boots < 0, so we set n.boots = 0', immediate. = TRUE)
    n.boots = 0
  }
  if (n.boots > 0) {
    extra.inputs = c('boot.ids', 'mlpmm.fits.boot')
    check2 = extra.inputs %in% ls(object)
    mess = paste('At least one of the following elements is missing in object:',
                 paste(extra.inputs, collapse = ', '))
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
  
  ########################
  ### original dataset ###
  ########################
  if (verbose) cat('Computing the predicted random effects on the original dataset...\n')
  # create matrix with predicted random effects computed on original dataset
  mlpmms = object$mlpmm.fits.orig
  
  ranef.orig = foreach(j = 1:p, .combine = 'cbind',
                       .packages = c('foreach', 'pencal', 'lcmm')) %dopar% {
    out = mlpmms[[j]]$predRE[,-1]
    rownames(out) = mlpmms[[j]]$predRE[,1]
    colnames(out) = paste(c('u0', 'u1'), j, sep = '_')
    if (randint.items[j]) {
      b = mlpmms[[j]]$predRE_Y[,-1]
      colnames(b) = paste('b0', j, colnames(b), sep = '_')
      out = cbind(out, b)
    }
    # return
    out
  }
  if (verbose) cat('...done\n')
  
  #######################
  ### start bootstrap ###
  #######################
  if (n.boots >= 1) {
    if (verbose) cat('Bootstrap procedure started\n')
    # retrieve bootstrap ids and fitted mlpmms
    boot.ids = object$boot.ids
    mlpmm.fits.boot = object$mlpmm.fits.boot
    
    # compute predicted ranefs for each bootstrap sample (training set)
    # and for the original dataset (validation set)
    
    # training set
    ranef.boot.train = foreach(t = 1:n.boots,
        .packages = c('foreach', 'pencal', 'lcmm')) %dopar% {
      ranef.train = foreach(j = 1:p, .combine = 'cbind') %do% {
        fit = mlpmm.fits.boot[[t]][[j]]
        out = fit$predRE[,-1]
        rownames(out) = fit$predRE[,1]
        colnames(out) = paste(c('u0', 'u1'), j, sep = '_')
        if (randint.items[j]) {
          b = fit$predRE_Y[,-1]
          colnames(b) = paste('b0', j, colnames(b), sep = '_')
          out = cbind(out, b)
        }
        # return
        out
      }
      # return
      ranef.train
    }
    
    # validation set
    # create fixed effects formulas for each group of items
    fixed.form = vector('list', p)
    for (i in 1:p) {
      f.left = paste(y.names[[i]], collapse = '+')
      fixed.form[[i]] = as.formula(paste(f.left, deparse(fixefs)))
    }
    
    # predict on the original dataset (called data) using the model 
    # fitted on the bootstrap sample
    ranef.boot.valid = foreach(t = 1:n.boots,
              .packages = c('foreach', 'pencal', 'nlme')) %dopar% {
      # for each response, derive predicted ranefs for that given bootstrap sample
      ranef.b = foreach(j = 1:p, .combine = 'cbind') %do% {
        fit = mlpmm.fits.boot[[t]][[j]]
        out = .ranef.mlpmm(fixed = fixed.form[[j]], random = ranefs, 
                          subject = 'numeric.id', 
                          mlpmm.fit = fit, newdata = df.sanitized)
        # return
        out
      }
      # fix the column names
      colnames(ranef.b) = names(ranef.boot.train[[t]])
      ranef.b
    }
    if (verbose) {
      cat('Bootstrap procedure finished\n')
      cat('Computation of step 2: finished :)\n')
    }
  }
  # close the cluster
  parallel::stopCluster(cl)
  
  out = list('call' = call, 'ranef.orig' = ranef.orig, 
             'n.boots' = n.boots)
  if (n.boots >= 1) {
    out[['boot.ids']] = boot.ids
    out[['ranef.boot.train']] = ranef.boot.train
    out[['ranef.boot.valid']] = ranef.boot.valid
  }
  class(out) = 'ranefs'
  return(out)
}
