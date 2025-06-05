#' Step 1 of PRC-MLPMM (estimation of the linear mixed models)
#'
#' This function performs the first step for the estimation
#' of the PRC-MLPMM model proposed in Signorelli et al. (2021)
#' 
#' @param y.names a list with the names of the
#' response variables which the MLPMMs have to be fitted to.
#' Each element in the list contains all the items used to 
#' reconstruct a latent biological process of interest
#' @param fixefs a fixed effects formula for the model, where the
#' time variable (specified also in \code{ranef.time}) is
#' included as first element and within the function 
#' \code{contrast()}. Examples: \code{~ contrast(age)}, 
#' \code{~ contrast(age) + group + treatment}
#' @param ranef.time a character with the name of the time variable 
#' for which to include a shared random slope
#' @param randint.items logical: should item-specific random intercepts
#' be included in the MLCMMs? Default is \code{TRUE}. It can also be a
#' vector, with different values for different elements of \code{y.names}
#' @param long.data a data frame with the longitudinal predictors,
#' comprehensive of a variable called \code{id} with the subject 
#' ids
#' @param surv.data a data frame with the survival data and (if 
#' relevant) additional baseline covariates. \code{surv.data} should at least
#' contain a subject id (called \code{id}), the time to event outcome  
#' (\code{time}), and binary event variable (\code{event})
#' @param t.from.base name of the variable containing time from 
#' baseline in \code{long.data}
#' @param n.boots number of bootstrap samples to be used in the
#' cluster bootstrap optimism correction procedure (CBOCP). If 0, no
#' bootstrapping is performed
#' @param n.cores number of cores to use to parallelize part of
#' the computations. If \code{ncores = 1} (default), 
#' no parallelization is done. Pro tip: you can use 
#' \code{parallel::detectCores()} to check how many 
#' cores are available on your computer
#' @param verbose if \code{TRUE} (default and recommended value), information
#' on the ongoing computations is printed in the console
#' @param seed random seed used for the bootstrap sampling. Default 
#' is \code{seed = 123}
#' @param maxiter maximum number of iterations to use when calling
#' the function \code{multlcmm}. Default is 100
#' @param conv a vector containing the three convergence criteria
#' (\code{convB}, \code{convL} and \code{convG}) to use when calling
#' the function \code{\link[lcmm]{multlcmm}}. Default is c(1e-3, 1e-3, 1e-3)
#' @param lcmm.warnings logical. If TRUE, a warning is printed every 
#' time the (strict) convergence criteria of the \code{multlcmm} function
#' are not met. Default is \code{FALSE}
#' 
#' @return A list containing the following objects:
#' \itemize{
#' \item \code{call.info}: a list containing the following function
#' call information: \code{call}, \code{y.names}, \code{fixefs},
#' \code{ranef.time}, \code{randint.items};
#' \item \code{mlpmm.fits.orig}: a list with the MLPMMs fitted on the
#' original dataset (it should comprise as many MLPMMs as the elements
#' of \code{y.names} are);
#' \item \code{df.sanitized}: a sanitized version of the supplied 
#' \code{long.data} dataframe, without the
#' longitudinal measurements that are taken after the event
#' or after censoring;
#' \item \code{n.boots}: number of bootstrap samples;
#' \item \code{boot.ids}: a list with the ids of bootstrapped subjects 
#' (when \code{n.boots > 0});
#' \item \code{mlpmm.fits.boot}: a list of lists, which contains the MLPMMs 
#' fitted on each bootstrapped datasets (when \code{n.boots > 0}).
#' }
#' 
#' @details This function is essentially a wrapper of the 
#' \code{\link[lcmm]{multlcmm}} that is meant to simplify
#' the estimation of several MLPMMs. In general, ensuring 
#' convergence of the algorithm implemented in \code{multlcmm}
#' is sometimes difficult, and it is hard to write a function that
#' can automatically solve all possible convergence problems. \code{fit_mplmms}
#' returns a warning when estimation did not converge for one or 
#' more MLPMMs. If this happens, try to change the convergence 
#' criteria in \code{conv} or the relevant \code{randint.items} value.
#' If doing this doesn't solve the problem, it is recommended to
#' re-estimate the specific MLPMMs for which estimation didn't converge
#' directly with \code{multlcmm}, trying to manually solve
#' the convergence issues
#' 
#' @import foreach doParallel stats
#' @importFrom lcmm multlcmm
#' @importFrom dplyr arrange
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
#' @seealso \code{\link{simulate_prcmlpmm_data}},
#' \code{\link{summarize_mlpmms}} (step 2),
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
#' # print MLPMM summary for marker 5 (all items involved in that MLPMM):
#' summary(step1, 'marker5_2')
#' }

fit_mlpmms = function(y.names, fixefs, ranef.time, 
                    randint.items = TRUE, long.data, 
                    surv.data, t.from.base, n.boots = 0, 
                    n.cores = 1, verbose = TRUE, seed = 123,
                    maxiter = 100, conv = rep(1e-3, 3),
                    lcmm.warnings = FALSE) {
  call = match.call()
  # load namespaces
  requireNamespace('lcmm')
  requireNamespace('foreach')
  requireNamespace('doParallel')
  # fix for 'no visible binding for global variable...' note
  id = i = numeric.id = b = NULL
  # setup and checks
  n.items = lengths(y.names)
  p = length(y.names)
  if (n.boots < 0) {
    warning('Input n.boots < 0, so we set n.boots = 0', immediate. = TRUE)
    n.boots = 0
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
  
  # checks on randint.items
  if (!is.logical(randint.items)) stop('randint.items has to be TRUE or FALSE')
  if (length(randint.items) == 1) {
    randint.items = rep(randint.items, length(y.names))
  }
  if (length(randint.items) != length(y.names)) {
    mess = paste('randint.items should either have length 1 or',
                  length(y.names))
    stop(mess)
  }
  # get ranef.time
  ranef.time = deparse(substitute(ranef.time))
  # get longitudinal predictors
  check1 = is.data.frame(long.data)
  if (!check1) stop('long.data should be a dataframe')
  check2 = c('id') %in% names(long.data)
  if (sum(check2) != 1) stop("long.data should contain an 'id' variable with subject ids")
  if (verbose) cat('Sorting long.data by subject id\n')
  long.data = dplyr::arrange(long.data, id)
  # get survival data
  check2 = is.data.frame(surv.data)
  if (!check2) stop('surv.data should be a dataframe')
  check3 = c('id', 'time', 'event') %in% names(surv.data)
  if (sum(check3) != 3) stop("surv.data should contain at least: 'id', 'time', 'event' variables")
  if (verbose) cat('Sorting surv.data by subject id\n')
  surv.data = dplyr::arrange(surv.data, id)
  # check that id values are the same in the two datasets!
  temp1 = as.character(unique(long.data$id))
  temp2 = as.character(unique(surv.data$id))
  check4 = identical(temp1, temp2)
  if (!check4) stop('id values are different in long.data and surv.data')  
  # remove all longitudinal measurements after t
  t.from.base = deparse(substitute(t.from.base))
  temp = prepare_longdata(df = long.data, subj.id = id, 
            t.from.base = t.from.base, survtime = surv.data$time, 
            verbose = verbose)
  df = temp$df.sanitized
  # add to df a numeric.id variable to simplify the bootstrapping
  df$numeric.id = as.numeric(as.factor(df$id))
  
  ########################
  ### original dataset ###
  ########################
  # set up environment for parallel computing
  cl = parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)
  .info_ncores(n.cores, verbose = verbose)
  
  # fit the MLPMMs
  if (verbose) cat('Estimating the MLPMMs on the original dataset...\n')
  ranef.formula = as.formula(paste('~ 1 +', ranef.time))
  # keep ranef.formula definition out of %dopar%
  
  fit.orig = foreach (i = 1:p, 
                .packages = c('foreach', 'pencal', 'lcmm')) %dopar% {
    # fit MLPMM
    ys = paste( paste(y.names[[i]], collapse = '+', sep = ''), sep = '')
    fixef.formula = as.formula(paste(ys, deparse(fixefs), sep = ''))

    mlpmm = try(multlcmm(fixed = fixef.formula, subject = 'id', 
             random = ranef.formula, randomY = randint.items[i],
             data = df, maxiter = maxiter, verbose = FALSE,
             convB = conv[1], convG = conv[2], convL = conv[3]))
    
    mess.part2 = paste('Try to increase maxiter, or (if this fails) to change the ',
                i, '-th element of randint.items. ',
                'See also point C of the Details section of ?multlcmm',
                sep = '')
    if (inherits(mlpmm, 'try-error')) {
      mess = paste('Estimation of the MLPMM for the ', i, '-th ',
             'latent biological process (on the original dataset) failed. ', 
             mess.part2,
             sep = '')
      stop(mess)
    }
    if (mlpmm$conv != 1 & lcmm.warnings) {
      mess = paste('Convergence not reached (on the original dataset) for the ', 
             i, '-th ', 'latent biological process. ', mess.part2, sep = '')
      warning(mess, immediate. = TRUE)
    }
    mlpmm
  }
  names(fit.orig) = y.names
  # close the cluster
  parallel::stopCluster(cl)
  if (verbose) cat('...done\n')
  
  #######################
  ### start bootstrap ###
  #######################
  if (n.boots >= 1) {
    if (verbose) cat('Bootstrap procedure started\n')
    # draw n bootstrap samples
    ids = unique(df$numeric.id)
    n = length(ids)
    set.seed(seed)
    boot.ids = foreach(i = 1:n.boots) %do% {
      sort(sample(ids, n, TRUE))
    }
    # set up environment for parallel computing
    cl = parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    
    # compute fitted LMMS for each bootstrap sample (in parallel)
    fit.boots = foreach(b = 1:n.boots, .export = 'draw_cluster_bootstrap',
    .packages = c('foreach', 'pencal', 'lcmm')) %dopar% {
      
      # bootstrap the longitudinal dataset for variable j
      boot.df = draw_cluster_bootstrap(df = df,
                        idvar = numeric.id, boot.ids = boot.ids[[b]])
      # important: overwrite id (that has repetitions) with new.id
      boot.df$id = boot.df$new.id
      
      # for each response, derive predicted ranefs for that given bootstrap sample
      this.fit = foreach(i = 1:p) %do% {
        # fit MLPMM
        ys = paste( paste(y.names[[i]], collapse = '+'))
        fixef.formula = as.formula(paste(ys, deparse(fixefs)))
        
        ranef.formula = as.formula(paste('~1+', ranef.time))
        
        mlpmm = try(multlcmm(fixed = fixef.formula, subject = 'id', 
                             random = ranef.formula, randomY = randint.items[i],
                             data = df, maxiter = maxiter, verbose = FALSE,
                             convB = conv[1], convG = conv[2], convL = conv[3]))
        
        if (inherits(mlpmm, 'try-error')) {
          mess = paste('Bootstrap sample', b,
                     ': estimation of the MLPMM for the ', i, '-th ',
                     'latent biological process failed. Try to increase maxiter', 
                     sep = '')
          warning(mess, immediate. = TRUE)
        }
        if (mlpmm$conv != 1 & lcmm.warnings) {
          mess = paste('Bootstrap sample', b,
                       ': convergence not reached for the ', i, '-th ',
                       'latent biological process. Try to increase maxiter', sep = '')
          warning(mess, immediate. = TRUE)
        }
        mlpmm
      }
      # all proteins done for each b
      names(this.fit) = y.names
      this.fit
    }
    # close the cluster
    parallel::stopCluster(cl)
    if (verbose) {
      cat('Bootstrap procedure finished\n')
      cat('Computation of step 1: finished :)\n')
    }
  }
  
  call.info = list('call' = call, 'y.names' = y.names,
                   'fixefs' = fixefs, 'ranefs' = ranef.formula,
                   'randint.items' = randint.items)
  out = list('call.info' = call.info, 'mlpmm.fits.orig' = fit.orig, 
             'df.sanitized' = df, 'n.boots' = n.boots)
  if (n.boots >= 1) {
    out[['boot.ids']] = boot.ids
    out[['mlpmm.fits.boot']] = fit.boots
  }
  class(out) = 'mlpmmfit'
  return(out)
}
