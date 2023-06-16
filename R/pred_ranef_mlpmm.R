.compute.sbj.ranef = function(sbj, df, id.variable,
  name.outcomes, link.mu, link.denom, name.fixefs,
  random, beta, contrs, sigma.u, sigma.b.vec,
  sigma.e.vec) {
  requireNamespace('magic')

  n.outcomes = length(name.outcomes)
  k = which(id.variable == sbj)
  df.sub = df[k, ]
  
  # create Yi
  Yi = c()
  for (i in 1:n.outcomes) {
    temp.Yi = as.numeric(df.sub[name.outcomes[i]][[1]])
    Yi = c(Yi, 
           (temp.Yi - link.mu[i])/link.denom[i] )
  }
  rm(temp.Yi)
  
  # create Xibeta and Zi
  fixef.formula = paste('~', paste(name.fixefs, collapse = '+'))
  Xi = model.matrix(as.formula(fixef.formula), df.sub)
  Zi.temp = model.matrix(random, df.sub)
  Zi = Zi.temp
  if (n.outcomes == 1) {
    Xibeta = Xi %*% beta
  }
  if (n.outcomes > 1) {
    beta.temp = beta + c(0, contrs[1])
    Xibeta = Xi %*% beta.temp
    for (i in 2:n.outcomes) {
      beta.temp = beta + c(0, contrs[i])
      Xibeta = rbind(Xibeta, Xi %*% beta.temp)
      # repeat Zi n.outcomes times
      Zi = rbind(Zi, Zi.temp)
    }
  }
  rm(Zi.temp)
  
  # first matrix in the computation
  el1 = t(Zi %*% sigma.u)
  el2 = kronecker(diag(sigma.b.vec), t(rep(1, length(k))))
  mat1 = rbind(el1, el2)
  
  # second matrix
  add1 = Zi %*% sigma.u %*% t(Zi)
  ntime = dim(Zi)[1]
  rtime = rep(length(k), n.outcomes)
  add2 = diag(rep(sigma.e.vec, rtime))
  add3 = do.call(magic::adiag, mapply(function(a,b)
    matrix(rep(a, b^2), ncol = b), sigma.b.vec, rtime, SIMPLIFY = FALSE))
  sigma.Yi = add1 + add2 + add3
  
  # third matrix
  res = Yi - Xibeta
  
  # computation BLUP
  blup = mat1 %*% solve(sigma.Yi) %*% res
  return(t(blup))
}

.ranef.mlpmm = function(fixed, random, subject = 'numeric.id', 
                      mlpmm.fit, newdata) {
  # fix for 'no visible binding for global variable...' note
  sbj = NULL
  
  # important: if there is more than 1 covariate in
  # fixed, make sure that fixed = ~ contrast(time) + x2 + x3
  # so that the time variable, with contrasts, is the very first one
  requireNamespace('lcmm')
  requireNamespace('foreach')
  n.outcomes = length(all.vars(fixed[[2]]))
  
  id.variable = newdata[, subject]
  numid = length(unique(id.variable))
  
  fit.names = names(mlpmm.fit$best)
  name.outcomes = all.vars(fixed[[2]])
  name.fixefs = all.vars(fixed[-2])
  
  # get positions
  pos.fixefs = which(fit.names %in% name.fixefs)
  pos.contrasts = which(substr(fit.names, 1, 8) == 'contrast')
  pos.sigma.u = which(substr(fit.names, 1, 6) == 'varcov')
  pos.sigma.e = which(substr(fit.names, 1, 7) == 'std.err')
  pos.sigma.b = which(substr(fit.names, 1, 11) == 'std.randomY')
  pos.link.mu = which(substr(fit.names, 1, 8) == 'Linear 1')
  pos.link.denom = which(substr(fit.names, 1, 8) == 'Linear 2')
  
  # retrieve variance matrices
  if (length(pos.sigma.u) != 2) {
    stop('code currently supports only 1 shared random intercept and slope')
  }
  if (length(pos.sigma.u) == 2) {
    sigma.u = matrix(c(1, mlpmm.fit$best['varcov 1'],
                       mlpmm.fit$best['varcov 1'], 
                       mlpmm.fit$best['varcov 2']), 2, 2)
  }
  sigma.e.vec = mlpmm.fit$best[pos.sigma.e]^2
  if (length(pos.sigma.b) > 0) sigma.b.vec = mlpmm.fit$best[pos.sigma.b]^2
  if (length(pos.sigma.b) == 0) sigma.b.vec = rep(0, n.outcomes)

  # get betas and links
  beta = c(0, mlpmm.fit$best[pos.fixefs])
  contrs = mlpmm.fit$best[pos.contrasts]
  last.contr = -sum(contrs)
  contrs = c(contrs, last.contr) # add the last element
  link.mu = mlpmm.fit$best[pos.link.mu]
  link.denom = mlpmm.fit$best[pos.link.denom]
  
  # compute estimates of random effects on new dataset
  id.variable = newdata[, subject]
  if (!is.numeric(id.variable)) id.variable = as.numeric(as.factor(id.variable))
  numid = length(unique(id.variable))
  pred.ranef = matrix(NA, nrow = numid, ncol = 2+n.outcomes)
  unique.sbjid = unique(id.variable)
    
  pred.ranef = foreach(sbj = unique.sbjid, .combine = rbind) %do%
    .compute.sbj.ranef(sbj = sbj, df = newdata, id.variable = id.variable,
                      name.outcomes = name.outcomes, link.mu = link.mu, 
                      link.denom = link.denom, name.fixefs = name.fixefs,
                      random = random, beta = beta, contrs = contrs, 
                      sigma.u = sigma.u, sigma.b.vec = sigma.b.vec, 
                      sigma.e.vec = sigma.e.vec)
  colnames(pred.ranef) = c('u0', 'u1',
                          paste('b', 1:n.outcomes, sep = ''))
  if (length(pos.sigma.b) == 0) {
    pred.ranef = pred.ranef[ , 1:2]
  }
  rownames(pred.ranef) = unique.sbjid
  pred.ranef = as.data.frame(pred.ranef)
  
  # define exports
  return(pred.ranef)
}
