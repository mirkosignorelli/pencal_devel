#' Generate survival data from a Weibull model
#'
#' This function implements the algorithm proposed by
#' Bender et al. (2005) to simulate survival times from 
#' a Weibull model
#' 
#' @param n sample size
#' @param lambda Weibull location parameter, positive
#' @param nu Weibull scale parameter, positive
#' @param X design matrix (n rows, p columns)
#' @param beta p-dimensional vector of regression coefficients
#' associated to X
#' @param seed random seed (defaults to 1)
#' 
#' @return A vector of survival times
#' 
#' @import stats
#' @export
#' 
#' @author Mirko Signorelli
#' @references 
#' Bender, R., Augustin, T., & Blettner, M. (2005). 
#' Generating survival times to simulate Cox proportional 
#' hazards models. Statistics in medicine, 24(11), 1713-1723.
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
#' set.seed(1)
#' n = 50
#' X = cbind(matrix(1, n, 1), 
#'    matrix(rnorm(n*9, sd = 0.7), n, 9))
#' beta = rnorm(10, sd = 0.7)
#' times = simulate_t_weibull(n = n, lambda = 1, nu = 2,
#'    X = X, beta = beta)
#' hist(times, 20)

simulate_t_weibull = function(n, lambda, nu, X, beta, seed = 1) {
  if (n <= 0) stop('n should be a positive integer')
  if (lambda <=0 | nu <= 0 | n <=0) stop('lambda, nu and n should be positive!')
  set.seed(seed)
  unif = runif(n, min = 0, max = 1)
  lin.pred = X %*% beta
  if (ncol(lin.pred) > 1) stop('error in computation of the linear predictor')
  temp = - log(unif) / (lambda * exp(lin.pred))
  time = temp^(1/nu)
  time = as.numeric(time)
  return(time)
}
