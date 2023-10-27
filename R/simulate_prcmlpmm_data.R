#' Simulate data that can be used to fit the PRC-LMM model
#'
#' This function allows to simulate a survival outcome
#' from longitudinal predictors following the PRC MLPMM model
#' presented in Signorelli et al. (2021). Specifically, the 
#' longitudinal predictors are simulated from multivariate 
#' latent process mixed models (MLPMMs), and 
#' the survival outcome from a Weibull model where the time
#' to event depends on the random effects from the MLPMMs.
#' 
#' @param n sample size
#' @param p number of longitudinal latent processes
#' @param p.relev number of latent processes that
#' are associated with the survival outcome (min: 1, max: p)
#' @param n.items number of items that are observed for each 
#' latent process of interest. It must be either a scalar, or
#' a vector of length \code{p}
#' @param type the type of relation between the longitudinal
#' outcomes and survival time. Two values can be used: 'u' 
#' refers to the PRC-MLPMM(U) model, and 'u+b' to the PRC-MLPMM(U+B)
#' model presented in Section 2.3 of Signorelli et al. (2021).
#' See the article for the mathematical details
#' @param t.values vector specifying the time points 
#' at which longitudinal measurements are collected
#' (NB: for simplicity, this function assumes a balanced 
#' designed; however, \code{pencal} is designed to work
#' both with balanced and with unbalanced designs!)
#' @param landmark the landmark time up until which all individuals survived.
#' Default is equal to \code{max(t.values)}
#' @param seed random seed (defaults to 1)
#' @param lambda Weibull location parameter, positive
#' @param nu Weibull scale parameter, positive
#' @param cens.range range for censoring times. By default, the minimum
#' of this range is equal to the \code{landmark} time
#' @param base.age.range range for age at baseline (set it
#' equal to c(0, 0) if you want all subjects to enter
#' the study at the same age)
#' @param tau.age the coefficient that multiplies baseline age
#' in the linear predictor (like in formulas (7) and (8) from  
#' Signorelli et al. (2021))
#' 
#' @return A list containing the following elements:
#' \itemize{
#' \item a dataframe \code{long.data} with data on the longitudinal 
#' predictors, comprehensive of a subject id (\code{id}),
#' baseline age (\code{base.age}), time from baseline
#' (\code{t.from.base}) and the longitudinal biomarkers;
#' \item a dataframe \code{surv.data} with the survival data: 
#' a subject id (\code{id}), baseline age (\code{baseline.age}),
#' the time to event outcome (\code{time}) and a binary vector
#' (\code{event}) that is 1 if the event
#' is observed, and 0 in case of right-censoring;
#' \item \code{perc.cens} the proportion of censored individuals 
#' in the simulated dataset.
#' }
#' 
#' @import stats
#' @importFrom Matrix bdiag
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
#' @examples
#' # generate example data
#' simdata = simulate_prcmlpmm_data(n = 40, p = 6,  
#'              p.relev = 3, n.items = c(3,4,2,5,4,2), 
#'              type = 'u+b', t.values = c(0, 0.5, 1, 2), 
#'              landmark = 2, seed = 19931101)
#' 
#' # names of the longitudinal outcomes:
#' names(simdata$long.data)
#' # markerx_y is the y-th item for latent process (LP) x
#' # we have 6 latent processes of interest, and for LP1 
#' # we measure 3 items, for LP2 4, for LP3 2 items, and so on
#' 
#' # visualize trajectories of marker1_1
#' if(requireNamespace("ptmixed")) {
#'   ptmixed::make.spaghetti(x = age, y = marker1_1, 
#'                  id = id, group = id,
#'                  data = simdata$long.data, 
#'                  legend.inset = - 1)
#'  }
#' # proportion of censored subjects
#' simdata$censoring.prop
#' # visualize KM estimate of survival
#' library(survival)
#' surv.obj = Surv(time = simdata$surv.data$time, 
#'                 event = simdata$surv.data$event)
#' kaplan <- survfit(surv.obj ~ 1,  
#'                  type="kaplan-meier")
#' plot(kaplan)

simulate_prcmlpmm_data = function(n = 100, p = 5, p.relev = 2,
              n.items = c(3, 2, 3, 4, 1), type = 'u',
              t.values = c(0, 0.5, 1, 2), landmark = max(t.values),
              seed = 1, lambda = 0.2, nu = 2, 
              cens.range = c(landmark, 10),
              base.age.range = c(3, 5), tau.age = 0.2) {
  if (n < 1) stop('n should be a positive integer')
  if (p < 1 | p.relev < 1) stop('p and p.relev should be a positive integer')
  if (p.relev > p) stop('p.relev should be <= p')
  if (lambda <=0 | nu <= 0 | n <=0) stop('lambda, nu and n should be positive!')
  # extra checks for the MLPMM:
  if (length(n.items) != p) {
    if (length(n.items) == 1) n.items = rep(n.items, p)
    else stop('n.items is neither of length 1 nor of length p')
  }
  check2 = all(n.items >=1, T)
  if (!check2) stop('all elements of n.items should be >=1')
  check3 = type %in% c('u', 'u+b')
  if (!check3) stop("vtype should be either 'u' or 'u+b'")
  
  requireNamespace('Matrix')
  
  set.seed(seed)
  n.total.items = sum(n.items)
  beta0 = runif(n.total.items, min = 3, max = 7)
  beta1 = .sim.eff.sizes(n = n.total.items, 
                        abs.range = c(1, 2))
  gamma = c(.sim.eff.sizes(n = p.relev,  
                          abs.range = c(0.5, 1),
                          seed = seed+1), 
            rep(0, p - p.relev))
  delta = c(.sim.eff.sizes(n = p.relev,  
                          abs.range = c(0.5, 1), 
                          seed = seed+2), 
                    rep(0, p - p.relev))
  if (type == 'u') xi = rep(0, n.total.items)
  if (type == 'u+b') {
    n.relev.items = sum(n.items[1:p.relev])
    xi = c(.sim.eff.sizes(n = n.relev.items,  
                          abs.range = c(0.2, 0.4), 
                          seed = seed+3), 
           rep(0, n.total.items - n.relev.items))
  }
  baseline.age = runif(n, base.age.range[1], base.age.range[2])

  # generate random effects
  sigma.b0 = diag(0.7, n.total.items)
  b0 = MASS::mvrnorm(n = n, mu = rep(0, n.total.items), Sigma = sigma.b0)
  sigma.u = matrix(c(1, 0.5, 0.5, 1), 2, 2)
  sigma = sigma.u
  for (i in 2:p) sigma = bdiag(sigma, sigma.u)
  us = MASS::mvrnorm(n = n, mu = rep(0, p*2), Sigma = sigma)
  u0 = us[, seq(1, 2*p-1, by = 2)]
  u1 = us[, seq(2, 2*p, by = 2)]
  
  m = length(t.values) # maximum theoretical number of repeated measurements
  id = rep(1:n, each = m)
  t.measures = rep(t.values, n) + rep(baseline.age, each = length(t.values))
  
  # create the longitudinal markers
  Y = matrix(NA, n*m, n.total.items)
  y.names = rep(NA, n.total.items)
  k = 1
  for (i in 1:p){ # latent process
    for (j in 1:n.items[i]) { # item for process i
      pos = j + ifelse(i ==1, 0, sum(n.items[1:(i-1)]))
      inter = rep(beta0[pos] + u0[,i] + b0[,pos], each = m)
      slope = rep(beta1[pos] + + u1[,i], each = m)
      mu = inter + slope*t.measures
      Y[,k] = rnorm(n = n*m, mean = mu, sd = 0.7)
      y.names[k] = paste('marker', i, '_', j, sep = '')
      k = k + 1
    }
  }
  
  long.data = cbind(id, rep(baseline.age, each = length(t.values)),
                    t.values, t.measures, Y)
  long.data = as.data.frame(long.data)
  names(long.data) = c('id', 'base.age', 't.from.base', 'age',
                       y.names)
  
  # simulate survival times
  X.temp = cbind(baseline.age, u0, u1, b0)
  true.t = simulate_t_weibull(n = n, lambda = lambda, nu = nu,
             beta = c(tau.age, gamma, delta, xi), 
             X = X.temp, seed = seed) + landmark
  censoring.times = runif(n = n, min = cens.range[1],
                        max = cens.range[2])
  event = ifelse(true.t > censoring.times, 0, 1)
  observed.t = ifelse(event == 1, true.t, censoring.times)
  surv.data = data.frame('id' = 1:n, baseline.age, 'time' = observed.t, event)
  censoring.prop = length(which(event == 0)) / length(event)
  out = list('long.data' = long.data, 'surv.data' = surv.data,
             'censoring.prop' = censoring.prop)
  return(out)
}

