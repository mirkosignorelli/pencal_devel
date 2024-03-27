#' Step 1 of PRC-LMM (estimation of the linear mixed models)
#'
#' This function performs the first step for the estimation
#' of the PRC-LMM model proposed in Signorelli et al. (2021)
#' 
#' @param y.names character vector with the names of the
#' response variables which the LMMs have to be fitted to
#' @param fixefs fixed effects formula for the model, example:
#' \code{~ time}
#' @param ranefs random effects formula for the model,
#' specified using the representation of random effect
#' structures of the \code{R} package \code{nlme}
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
#' @param max.ymissing maximum proportion of subjects allowed to not have any
#' measurement of a longitudinal response variable. Default is 0.2
#' @param verbose if \code{TRUE} (default and recommended value), information
#' on the ongoing computations is printed in the console
#' @param seed random seed used for the bootstrap sampling. Default 
#' is \code{seed = 123}
#' @param control a list of control values to be passed to \code{lme} when fitting the
#' linear mixed models. By default, we set \code{opt = 'optim', niterEM = 500, maxIter = 500}. 
#' See \code{?nlme::lmeControl} for all possible arguments and values
#' 
#' @return A list containing the following objects:
#' \itemize{
#' \item \code{call.info}: a list containing the following function
#' call information: \code{call}, \code{y.names}, \code{fixefs},
#' \code{ranefs};
#' \item \code{lmm.fits.orig}: a list with the LMMs fitted on the
#' original dataset (it should comprise as many LMMs as the elements
#' of \code{y.names} are);
#' \item \code{df.sanitized}: a sanitized version of the supplied 
#' \code{long.data} dataframe, without the
#' longitudinal measurements that are taken after the event
#' or after censoring;
#' \item \code{n.boots}: number of bootstrap samples;
#' \item \code{boot.ids}: a list with the ids of bootstrapped subjects 
#' (when \code{n.boots > 0});
#' \item \code{lmms.fits.boot}: a list of lists, which contains the LMMs fitted 
#' on each bootstrapped datasets (when \code{n.boots > 0}).
#' }
#' 
#' @import nlme foreach doParallel stats
#' @importFrom dplyr arrange
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
#' @seealso \code{\link{simulate_prclmm_data}},
#' \code{\link{summarize_lmms}} (step 2),
#' \code{\link{fit_prclmm}} (step 3),
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
#' # estimated betas and variances for the 3rd marker:
#' summary(step1, 'marker3', 'betas')
#' summary(step1, 'marker3', 'variances')
#' # usual T table:
#' summary(step1, 'marker3', 'tTable')

fit_lmms = function(y.names, fixefs, ranefs, long.data, 
                    surv.data, t.from.base, n.boots = 0, 
                    n.cores = 1, max.ymissing = 0.2, 
                    verbose = TRUE, seed = 123,
                    control = list(opt = 'optim', niterEM = 500, maxIter = 500)) {
  call = match.call()
  # load namespaces
  requireNamespace('nlme')
  requireNamespace('foreach')
  requireNamespace('doParallel')
  # fix for 'no visible binding for global variable...' note
  id = i = numeric.id = b = NULL
  # setup and checks
  p = length(y.names)
  if (n.boots < 0) {
    warning('Input n.boots < 0, so we set n.boots = 0', immediate. = TRUE)
    n.boots = 0
  }
  if (max.ymissing < 0 | max.ymissing >= 1) {
    stop('max.ymissing should be in [0, 1)')
  }
  if (max.ymissing > 0.4) {
    warning(paste('You have set ymissing = ', max.ymissing,
               '. Are you sure about that?', sep = ''), immediate. = TRUE)
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
  
  # set up environment for parallel computing
  cl = parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)
  .info_ncores(n.cores, verbose = verbose)

  ########################
  ### original dataset ###
  ########################
  # fit the LMMs
  if (verbose) cat('Estimating the LMMs on the original dataset...\n')
  fit.orig = foreach(i = 1:p,
    .packages = c('foreach', 'pencal', 'nlme')) %dopar% {
    # check if NAs on response
    n1 = length(unique(df$numeric.id))
    nas = which(is.na(df[ , y.names[i]]))
    df.sub = df
    if (length(nas) > 0) df.sub = df[-nas, ]
    n2 = length(unique(df.sub$numeric.id))
    if (n1 != n2) {
      if (max.ymissing == 0) {
        mess = paste('There is at least one subject without any information available for variable ',
                     y.names[i], '. Consider increasing the value of max.ymissing', sep = '')
        stop(mess)
      }
      prop.na = 1 - n2/n1
      if (max.ymissing < prop.na) {
        mess = paste('The proportion of subject with missing ', y.names[i],
                     ' is above ', max.ymissing, '. Either increase max.ymissing',
                     ' or remove ', y.names[i], ' from the longitudinal modelling', sep = '')
        stop(mess)
      }
    }
    # fit LMM
    fixef.formula = as.formula(paste(y.names[i], deparse(fixefs)))
    ranef.formula = as.formula(deparse(ranefs))
    lmm = try( nlme::lme(fixed = fixef.formula, 
                    random = ranef.formula, 
                    data = df, na.action = na.exclude,
                    keep.data = FALSE, control = control),
                silent = TRUE)
    if (inherits(lmm, 'try-error')) { 
      # retry with larger max number of iterations (larger than lmeControl's controls given above)
      lmm = try( nlme::lme(fixed = fixef.formula, 
                     random = ranef.formula, 
                     data = df, na.action = na.exclude,
                     keep.data = FALSE,
                     control = list(control$opt, niterEM = 1e3,
                                    maxIter = 1e3, msMaxIter = 1e3, msMaxEval = 1e3)),
                 silent = TRUE)
    }
    if (inherits(lmm, 'try-error')) {
      mess = paste('The desired model for response', y.names[i],
                   'did not converge. A simpler model comprising only a random
                   intercept was fitted for', y.names[i])
      warning(mess, immediate. = TRUE)
      lmm = try( nlme::lme(fixed = fixef.formula, 
                           random = ~ 1 | id, 
                           data = df, na.action = na.exclude,
                           keep.data = FALSE,
                           control = control),
                 silent = TRUE)
    }
    if (inherits(lmm, 'try-error')) {
      stop(paste('no model could be fitted for response', y.names[i]))
    }
    lmm
  }
  names(fit.orig) = y.names
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
    
    # compute fitted LMMS for each bootstrap sample (in parallel)
    fit.boots = foreach(b = 1:n.boots,
      .packages = c('foreach', 'pencal', 'nlme')) %dopar% {
      
      # bootstrap the longitudinal dataset for variable j
      boot.df = draw_cluster_bootstrap(df = df,
                        idvar = numeric.id, boot.ids = boot.ids[[b]])
      # important: overwrite id (that has repetitions) with new.id
      boot.df$id = boot.df$new.id
      
      # for each response, derive predicted ranefs for that given bootstrap sample
      this.fit = foreach(i = 1:p) %do% {
        # fit LMM
        fixef.formula = as.formula(paste(y.names[i], deparse(fixefs)))
        ranef.formula = as.formula(deparse(ranefs))
        lmm = try( nlme::lme(fixed = fixef.formula, 
                              random = ranef.formula, 
                              data = boot.df, na.action = na.exclude,
                              keep.data = FALSE, control = control),
                    silent = TRUE)
        if (inherits(lmm, 'try-error')) {
          # retry with larger max number of iterations (larger than lmeControl's defaults)
          lmm = try( nlme::lme(fixed = fixef.formula, 
                               random = ranef.formula, 
                               data = boot.df, na.action = na.exclude,
                               keep.data = FALSE,
                               control = list(control$opt, niterEM = 1e3,
                                              maxIter = 1e3, msMaxIter = 1e3, 
                                              msMaxEval = 1e3),
                     silent = TRUE))
        }
        if (inherits(lmm, 'try-error')) {
          lmm = try( nlme::lme(fixed = fixef.formula, 
                               random = ~ 1 | id, 
                               data = boot.df, na.action = na.exclude,
                               keep.data = FALSE,
                               control = control),
                     silent = TRUE)
        }
        if (inherits(lmm, 'try-error')) {
          warning(paste('Bootstrap sample', b,
            ': the LMM could not be fitted for response', y.names[i]), 
            immediate. = TRUE)
        }
        lmm
      }
      # all Ys done for each bootstrap sample
      names(this.fit) = y.names
      this.fit
    }
    if (verbose) {
      cat('Bootstrap procedure finished\n')
      cat('Computation of step 1: finished :)\n')
    }
  }
  # close the cluster
  parallel::stopCluster(cl)
  
  # define exports
  call.info = list('call' = call, 'y.names' = y.names,
                   'fixefs' = fixefs, 'ranefs' = ranefs)
  out = list('call.info' = call.info, 'lmm.fits.orig' = fit.orig, 
             'df.sanitized' = df, 'n.boots' = n.boots)
  if (n.boots >= 1) {
    out[['boot.ids']] = boot.ids
    out[['lmm.fits.boot']] = fit.boots
  }
  class(out) = 'lmmfit'
  return(out)
}
