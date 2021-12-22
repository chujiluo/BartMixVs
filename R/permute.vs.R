permute.vs = function(x.train, 
                      y.train, 
                      probit=F, 
                      npermute=100L,                      ## number of permutations 
                      nreps=10L,                          ## number of replicates 
                      alpha=0.05,                         ## local threshold
                      true.idx=NULL,
                      plot=T, 
                      n.var.plot = Inf,
                      xinfo=matrix(0.0,0,0), 
                      numcut=100L,
                      usequants=FALSE, 
                      cont=FALSE, 
                      rm.const=TRUE, 
                      k=2.0,
                      power=2.0, 
                      base=0.95,
                      split.prob="polynomial",
                      ntree=20L, 
                      ndpost=1000, 
                      nskip=1000,
                      keepevery=1L, 
                      printevery=100L,
                      verbose=FALSE) {
  
  res = list()
  varcounts = list() ## collect all the varcount matrices
  
  #------------------------------
  # timer starts
  start = Sys.time()
  
  #-----------------------------------------------------------
  # data
  categorical.idx = which(sapply(x.train, function(s) {is.factor(s)}))
  categorical.names = names(categorical.idx)
  
  #-----------------------------------------------------------
  # get avg/median variable importance from the original data
  cat("original data set...")
  
  avg.vip.mtx = matrix(NA, nrow = nreps, ncol = ncol(x.train))
  median.mi.mtx = matrix(NA, nrow = nreps, ncol = ncol(x.train))
  if (length(categorical.idx) > 0)
    avg.within.type.vip.mtx = matrix(NA, nrow = nreps, ncol = ncol(x.train))
  
  cnt = 0
  while (cnt < nreps) {
    if (probit) {
      bart = pbart(x.train = x.train, y.train = y.train, sparse = FALSE,
                   xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, rm.const = rm.const,
                   k = k, power = power, base = base, split.prob = split.prob,
                   ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, verbose = verbose) 
    } else {
      bart = wbart(x.train = x.train, y.train = y.train, sparse = FALSE,
                   xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, rm.const = rm.const,
                   k = k, power = power, base = base, split.prob = split.prob,
                   ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, verbose = verbose) 
    }
    
    cnt = cnt + 1
    varcounts[[cnt]] = bart$varcount
    avg.vip.mtx[cnt, ] = bart$vip
    median.mi.mtx[cnt, ] = bart$mi
    if (length(categorical.idx) > 0)
      avg.within.type.vip.mtx[cnt, ] = bart$within.type.vip
  }
    
  avg.vip = colMeans(avg.vip.mtx)
  names(avg.vip) = colnames(x.train)
  avg.vip = sort(avg.vip, decreasing = T)
  
  median.mi = apply(median.mi.mtx, 2, median)
  names(median.mi) = colnames(x.train)
  median.mi = sort(median.mi, decreasing = T)
  
  if (length(categorical.idx) > 0) {
    avg.within.type.vip = colMeans(avg.within.type.vip.mtx)
    names(avg.within.type.vip) = colnames(x.train)
    avg.within.type.vip = sort(avg.within.type.vip, decreasing = T)
  }

  cat("complete! \n")
  
  
  #-----------------------------------------------------------
  # build null permutation
  cat("null data sets...")
  
  ## set up permute matrix
  permute.vips = matrix(NA, nrow = npermute, ncol = ncol(x.train))
  permute.mis = matrix(NA, nrow = npermute, ncol = ncol(x.train))
  if (length(categorical.idx) > 0)
    permute.within.type.vips = matrix(NA, nrow = npermute, ncol = ncol(x.train))
  
  cnt = 0
  while (cnt < npermute) {
    y.permuted = sample(y.train, replace = F)
    
    if (probit) {
      bart = pbart(x.train = x.train, y.train = y.permuted, sparse = FALSE,
                   xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, rm.const = rm.const,
                   k = k, power = power, base = base, split.prob = split.prob,
                   ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, verbose = verbose) 
    } else {
      bart = wbart(x.train = x.train, y.train = y.permuted, sparse = FALSE,
                   xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, rm.const = rm.const,
                   k = k, power = power, base = base, split.prob = split.prob,
                   ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, verbose = verbose) 
    }
    
    cnt = cnt + 1
    
    varcounts[[nreps + cnt]] = bart$varcount
    permute.vips[cnt, ] = bart$vip
    permute.mis[cnt, ] = bart$mi
    if (length(categorical.idx) > 0)
      permute.within.type.vips[cnt, ] = bart$within.type.vip
  }
  
  cat("complete! \n")
  
  
  #-----------------------------------------------------------
  # sort permute mat and return results
  colnames(permute.vips) = colnames(x.train)
  permute.vips = permute.vips[, names(avg.vip)]
  
  colnames(permute.mis) = colnames(x.train)
  permute.mis = permute.mis[, names(median.mi)]
  
  if (length(categorical.idx) > 0) {
    colnames(permute.within.type.vips) = colnames(x.train)
    permute.within.type.vips = permute.within.type.vips[, names(avg.within.type.vip)]
  }
  
  #-----------------------------------------------------------
  # use local cutoff & returns
  vip.pointwise.cutoffs = apply(permute.vips, 2, quantile, probs = 1 - alpha)
  vip.imp.names = names(avg.vip[(avg.vip > vip.pointwise.cutoffs) & (avg.vip > 0)])
  vip.imp.cols = sapply(1:length(vip.imp.names), function(x) {which(vip.imp.names[x] == colnames(x.train))})
  res$vip.imp.cols = vip.imp.cols
  res$vip.imp.names = vip.imp.names
  res$avg.vip = avg.vip
  res$avg.vip.mtx = avg.vip.mtx
  res$permute.vips = permute.vips
  
  mi.pointwise.cutoffs = apply(permute.mis, 2, quantile, probs = 1 - alpha)
  mi.imp.names = names(median.mi[(median.mi > mi.pointwise.cutoffs) & (median.mi > 0)])
  mi.imp.cols = sapply(1:length(mi.imp.names), function(x) {which(mi.imp.names[x] == colnames(x.train))})
  res$mi.imp.cols = mi.imp.cols
  res$mi.imp.names = mi.imp.names
  res$median.mi = median.mi
  res$median.mi.mtx = median.mi.mtx
  res$permute.mis = permute.mis
  
  if (length(categorical.idx) > 0) {
    within.type.vip.pointwise.cutoffs = apply(permute.within.type.vips, 2, quantile, probs = 1 - alpha)
    within.type.vip.imp.names = names(avg.within.type.vip[(avg.within.type.vip > within.type.vip.pointwise.cutoffs) & (avg.within.type.vip > 0)])
    within.type.vip.imp.cols = sapply(1:length(within.type.vip.imp.names), function(x) {which(within.type.vip.imp.names[x] == colnames(x.train))})
    
    res$within.type.vip.imp.cols = within.type.vip.imp.cols
    res$within.type.vip.imp.names = within.type.vip.imp.names
    res$avg.within.type.vip = avg.within.type.vip
    res$avg.within.type.vip.mtx = avg.within.type.vip.mtx
    res$permute.within.type.vips = permute.within.type.vips
  }
  
  res$varcounts = varcounts
  
  #-----------------------------------------------------------
  # score results
  if(length(true.idx) > 0) {
    
    true.len = length(true.idx)
    res$true.idx = true.idx
    
    ## vip
    tp = length(which(vip.imp.cols %in% true.idx))
    positive.len = length(vip.imp.cols)
    
    res$vip.precision = (tp * 1.0) / (positive.len * 1.0)
    res$vip.recall = (tp * 1.0) / (true.len * 1.0)
    res$vip.f1 = 2 * res$vip.precision * res$vip.recall / (res$vip.precision + res$vip.recall)
    
    ## mi
    tp = length(which(mi.imp.cols %in% true.idx))
    positive.len = length(mi.imp.cols)
    
    res$mi.precision = (tp * 1.0) / (positive.len * 1.0)
    res$mi.recall = (tp * 1.0) / (true.len * 1.0)
    res$mi.f1 = 2 * res$mi.precision * res$mi.recall / (res$mi.precision + res$mi.recall)
    
    if (length(categorical.idx) > 0) {
      ## within-type vip
      tp = length(which(within.type.vip.imp.cols %in% true.idx))
      positive.len = length(within.type.vip.imp.cols)
      
      res$wt.vip.precision = (tp * 1.0) / (positive.len * 1.0)
      res$wt.vip.recall = (tp * 1.0) / (true.len * 1.0)
      res$wt.vip.f1 = 2 * res$wt.vip.precision * res$wt.vip.recall / (res$wt.vip.precision + res$wt.vip.recall)
    }
  }
  
  #-----------------------------------------------------------
  if (plot) {
    
    if ((n.var.plot == Inf) | (n.var.plot > ncol(x.train))){
      n.var.plot = ncol(x.train)
    }
    
    ## vip
    non.zero.idx = which(avg.vip > 0)[1:min(n.var.plot, length(which(avg.vip > 0)))]
    plot.n = length(non.zero.idx)
    if(length(non.zero.idx) < length(avg.vip)) 
      warning(paste(length(which(avg.vip == 0)), "predictors with inclusion proportions of 0 omitted from plots."))
    max.cut = max(apply(permute.vips, 2, quantile, probs = 1 - alpha, na.rm = TRUE))
    
    plot(1:plot.n, avg.vip[non.zero.idx], 
         type = "n", xlab = "Predictors", xaxt = "n", ylim = c(0, max(max(avg.vip), max.cut * 1.1)),
         main = "Permutation-Based Variable Selection", ylab = "BART VIP")
    axis(1, at = 1:plot.n, labels = names(avg.vip[non.zero.idx]), las = 2)
    for (j in non.zero.idx){
      points(j, avg.vip[j], 
             pch = ifelse(avg.vip[j] <= quantile(permute.vips[, j], 1 - alpha), 1, 16),
             col = ifelse(names(avg.vip[j]) %in% categorical.names, 'green', 'red'))
    }
    sapply(non.zero.idx, function(s) {segments(s, 0, x1 = s, quantile(permute.vips[, s], 1 - alpha), col = "grey")})
    
    ## mi
    non.zero.idx = which(median.mi > 0)[1:min(n.var.plot, length(which(median.mi > 0)))]
    plot.n = length(non.zero.idx)
    if(length(non.zero.idx) < length(median.mi)) 
      warning(paste(length(which(median.mi == 0)), "predictors with Metropolis importance of 0 omitted from plots."))
    max.cut = max(apply(permute.mis, 2, quantile, probs = 1 - alpha, na.rm = TRUE))
    
    plot(1:plot.n, median.mi[non.zero.idx], 
         type = "n", xlab = "Predictors", xaxt = "n", ylim = c(0, max(max(median.mi), max.cut * 1.1)),
         main = "Permutation-Based Variable Selection", ylab = "BART MI")
    axis(1, at = 1:plot.n, labels = names(median.mi[non.zero.idx]), las = 2)
    for (j in non.zero.idx){
      points(j, median.mi[j], 
             pch = ifelse(median.mi[j] <= quantile(permute.mis[, j], 1 - alpha), 1, 16),
             col = ifelse(names(median.mi[j]) %in% categorical.names, 'green', 'red'))
    }
    sapply(non.zero.idx, function(s) {segments(s, 0, x1 = s, quantile(permute.mis[, s], 1 - alpha), col = "grey")})
    
    if (length(categorical.idx) > 0) {
      ## within-type vip
      non.zero.idx = which(avg.within.type.vip > 0)[1:min(n.var.plot, length(which(avg.within.type.vip > 0)))]
      plot.n = length(non.zero.idx)
      if(length(non.zero.idx) < length(avg.within.type.vip)) 
        warning(paste(length(which(avg.within.type.vip == 0)), 
                      "predictors with inclusion proportions (within-type) of 0 omitted from plots."))
      max.cut = max(apply(permute.within.type.vips, 2, quantile, probs = 1 - alpha, na.rm = TRUE))
      
      plot(1:plot.n, avg.within.type.vip[non.zero.idx], 
           type = "n", xlab = "Predictors", xaxt = "n", ylim = c(0, max(max(avg.within.type.vip), max.cut * 1.1)),
           main = "Permutation-Based Variable Selection", ylab = "BART Within-Type VIP")
      axis(1, at = 1:plot.n, labels = names(avg.within.type.vip[non.zero.idx]), las = 2)
      for (j in non.zero.idx){
        points(j, avg.within.type.vip[j], 
               pch = ifelse(avg.within.type.vip[j] <= quantile(permute.within.type.vips[, j], 1 - alpha), 1, 16),
               col = ifelse(names(avg.within.type.vip[j]) %in% categorical.names, 'green', 'red'))
      }
      sapply(non.zero.idx, function(s) {segments(s, 0, x1 = s, quantile(permute.within.type.vips[, s], 1 - alpha), col = "grey")})
    }
  }
  
  #------------------------------
  # timer ends
  end = Sys.time()
  cat("Elapsed", end-start, '\n')
  
  return(res)
}
               
