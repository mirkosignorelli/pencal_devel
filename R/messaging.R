.check_ncores = function(requested, avail, verbose = TRUE) {
  if (verbose) {
    diff = avail - requested
    mess0 = paste('You requested', requested, 'cores for this computation.')
    mess1 = paste(mess0, 'It seems that your computer actually has',
                  avail, 'cores available.',
                  'Consider increasing n.cores accordingly to speed computations up! =)')
    if (diff > 0) message(mess1)
    mess2 = paste(mess0, 'However, it seems that your computer only has',
                  avail, 'cores available.',
                  'Therefore, computations will be performed using only', 
                  avail, 'cores. =(')
    if (diff < 0)  message(mess2) 
  }
}

.info_ncores = function(n.cores, verbose = TRUE) {
  mess = ifelse(n.cores >= 2, 'in parallel',
                'without parallelization')
  if (verbose) cat(paste('This computation will be run', mess, '\n'))
}
