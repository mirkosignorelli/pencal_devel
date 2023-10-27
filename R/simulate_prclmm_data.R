#' Simulate data that can be used to fit the PRC-LMM model
#'
#' This function allows to simulate a survival outcome
#' from longitudinal predictors following the PRC LMM model
#' presented in Signorelli et al. (2021). Specifically, the longitudinal
#' predictors are simulated from linear mixed models (LMMs), and 
#' the survival outcome from a Weibull model where the time
#' to event depends linearly on the baseline age and on the 
#' random effects from the LMMs.
#' 
#' @param n sample size
#' @param p number of longitudinal outcomes
#' @param p.relev number of longitudinal outcomes that
#' are associated with the survival outcome (min: 1, max: p)
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
#' in the linear predictor (like in formula (6) from Signorelli 
#' et al. (2021))
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
#' in the simulated dataset;
#' \item \code{theta.true} a list containing the true parameter values
#' used to simulate data from the mixed model (beta0 and beta1) and
#' from the Weibull model (tau.age, gamma, delta)
#' }
#' 
#' @import stats
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
#' simdata = simulate_prclmm_data(n = 20, p = 10, p.relev = 4,
#'                t.values = c(0, 0.5, 1, 2), landmark = 2, 
#'                seed = 19931101)
#' # view the longitudinal markers:
#' if(requireNamespace("ptmixed")) {
#'   ptmixed::make.spaghetti(x = age, y = marker1, 
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
#'                   type="kaplan-meier")
#' plot(kaplan)

simulate_prclmm_data = function(n = 100, p = 10, p.relev = 4,
              t.values = c(0, 0.5, 1, 2), landmark = max(t.values),
              seed = 1, lambda = 0.2, nu = 2, cens.range = c(landmark, 10),
              base.age.range = c(3, 5), tau.age = 0.2) {
  if (n < 1) stop('n should be a positive integer')
  if (p < 1 | p.relev < 1) stop('p and p.relev should be a positive integer')
  if (p.relev > p) stop('p.relev should be <= p')
  if (lambda <=0 | nu <= 0 | n <=0) stop('lambda, nu and n should be positive!')

  set.seed(seed)
  beta0 = runif(p, min = 3, max = 7)
  beta1 = .sim.eff.sizes(n = p, abs.range = c(1, 2))
  gamma = c(.sim.eff.sizes(n = p.relev,  
                          abs.range = c(0.5, 1),
                          seed = seed + 1), 
            rep(0, p - p.relev))
  delta = c(.sim.eff.sizes(n = p.relev,  
                          abs.range = c(0.5, 1), 
                          seed = seed + 2), 
                    rep(0, p - p.relev))
  baseline.age = runif(n, base.age.range[1], base.age.range[2])

  Sigma = diag(1, p)
  rand.int = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  rand.slope = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  m = length(t.values) # maximum theoretical number of repeated measurements
  id = rep(1:n, each = m)
  t.measures = rep(t.values, n) + rep(baseline.age, each = length(t.values))
  
  Y = matrix(NA, n*m, p)
  
  for (i in 1:p){
    inter = rep(beta0[i] + rand.int[,i], each = m)
    slope = rep(beta1[i] + rand.slope[,i], each = m)
    mu = inter + slope*t.measures
    Y[,i] = rnorm(n = n*m, mean = mu, sd = 0.7)
  }
  
  long.data = cbind(id, rep(baseline.age, each = length(t.values)),
                    t.values, t.measures, Y)
  long.data = as.data.frame(long.data)
  names(long.data) = c('id', 'base.age', 't.from.base', 'age',
                       paste('marker', 1:p, sep = ''))
  
  # simulate survival times
  X.temp = cbind(baseline.age, rand.int, rand.slope)
  true.t = simulate_t_weibull(n = n, lambda = lambda, nu = nu,
             beta = c(tau.age, gamma, delta), 
             X = X.temp, seed = seed) + landmark
  censoring.times = runif(n = n, min = cens.range[1],
                        max = cens.range[2])
  event = ifelse(true.t > censoring.times, 0, 1)
  observed.t = ifelse(event == 1, true.t, censoring.times)
  surv.data = data.frame('id' = 1:n, baseline.age, 'time' = observed.t, event)
  censoring.prop = length(which(event == 0)) / length(event)
  theta.true = list('beta0' = beta0, 'beta1' = beta1,
               'tau.age' = tau.age, 'gamma' = gamma, 'delta' = delta)
  out = list('long.data' = long.data, 'surv.data' = surv.data,
             'censoring.prop' = censoring.prop, 'theta.true' = theta.true)
  return(out)
}

.sim.eff.sizes = function(n, abs.range, seed = 1) {
  # simulates from a U(a,b) and adds a sign (p = 0.5)
  set.seed(seed)
  abs.val = stats::runif(n, min = abs.range[1], max = abs.range[2])
  signs = 2*stats::rbinom(n, 1, prob = 0.5) - 1
  out = signs*abs.val
  return(out)
}
