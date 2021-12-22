median.inclusion.vs = function(x.train, 
                               y.train, 
                               probit=FALSE,             ## pbart or wbart    
                               vip.selection=TRUE,       ## if do VS using BART VIP with threshold = 1/p
                               true.idx=NULL, 
                               plot=F, 
                               num.var.plot=Inf,
                               theta=0,                  ## DART parameters
                               omega=1,                      
                               a=0.5, 
                               b=1, 
                               augment=FALSE, 
                               rho=NULL,
                               xinfo=matrix(0.0,0,0), 
                               numcut=100L,
                               usequants=FALSE, 
                               cont=FALSE, 
                               rm.const=TRUE, 
                               power=2.0, 
                               base=.95,                   ## depth parameters
                               split.prob="polynomial",
                               k=2.0,                      ## leaf parameters
                               ntree=20L, 
                               ndpost=1000L, 
                               nskip=1000L, 
                               keepevery=1L, 
                               printevery=100L,
                               verbose=FALSE) {
  
  res = list()
  
  #------------------------------
  # timer starts
  start = Sys.time()
  
  #------------------------------------------
  # fit bart and/or dart
  if (probit) {
    dart = pbart(x.train = x.train, y.train = y.train,
                 sparse = TRUE, theta = theta, omega = omega, a = a, b = b, augment = augment, rho = rho,
                 xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, rm.const = rm.const, 
                 power = power, base = base, split.prob = split.prob, k = k, 
                 ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, 
                 printevery = printevery, verbose = verbose)
    
    if (vip.selection) {
      bart = pbart(x.train = x.train, y.train = y.train, sparse = FALSE, 
                   xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, rm.const = rm.const, 
                   power = power, base = base, split.prob = split.prob, k = k, 
                   ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, 
                   printevery = printevery, verbose = verbose)
    }
  } else {
    dart = wbart(x.train = x.train, y.train = y.train,
                 sparse = TRUE, theta = theta, omega = omega, a = a, b = b, augment = augment, rho = rho,
                 xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, rm.const = rm.const, 
                 power = power, base = base, split.prob = split.prob, k = k, 
                 ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, 
                 printevery = printevery, verbose = verbose)
    
    if (vip.selection) {
      bart = wbart(x.train = x.train, y.train = y.train, sparse = FALSE, 
                   xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, rm.const = rm.const, 
                   power = power, base = base, split.prob = split.prob, k = k, 
                   ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, 
                   printevery = printevery, verbose = verbose)
    }
  }
  
  #------------------------------------------
  # get importance scores
  
  ## DART marginal posterior inclusion probability, threshold = 0.5
  dart.pvip = dart$pvip
  dart.pvip = sort(dart.pvip, decreasing = T)
  dart.pvip.imp.names = names(dart.pvip[dart.pvip > 0.5])
  dart.pvip.imp.cols = sapply(1:length(dart.pvip.imp.names),
                              function(x) {which(dart.pvip.imp.names[x] == colnames(x.train))})
  
  res$dart.pvip = dart.pvip
  res$dart.pvip.imp.names = dart.pvip.imp.names
  res$dart.pvip.imp.cols = dart.pvip.imp.cols
  
  rm(dart)
  
  if (vip.selection) {
    ## BART variable inclusion proportions, threshold = 1/p
    bart.vip = bart$vip
    bart.vip = sort(bart.vip, decreasing = T)
    bart.vip.imp.names = names(bart.vip[bart.vip > (1/dim(x.train)[2])])
    bart.vip.imp.cols = sapply(1:length(bart.vip.imp.names),
                               function(x) {which(bart.vip.imp.names[x] == colnames(x.train))})
    
    res$bart.vip = bart.vip
    res$bart.vip.imp.names = bart.vip.imp.names
    res$bart.vip.imp.cols = bart.vip.imp.cols
    
    rm(bart)
  }
  
  
  #------------------------------------------
  # score the results
  if (length(true.idx) > 0) {
    true.len = length(true.idx)
    
    dart.tp = length(which(dart.pvip.imp.cols %in% true.idx))
    dart.positive.len = length(dart.pvip.imp.cols)
    
    dart.prec = (dart.tp * 1.0) / (dart.positive.len * 1.0)
    dart.rec = (dart.tp * 1.0) / (true.len * 1.0)
    dart.f1 = 2 * dart.prec * dart.rec / (dart.prec + dart.rec)
    
    res$dart.precision = dart.prec
    res$dart.recall = dart.rec
    res$dart.f1 = dart.f1
    
    if(vip.selection) {
      bart.tp = length(which(bart.vip.imp.cols %in% true.idx))
      bart.positive.len = length(bart.vip.imp.cols)
      
      bart.prec = (bart.tp * 1.0) / (bart.positive.len * 1.0)
      bart.rec = (bart.tp * 1.0) / (true.len * 1.0)
      bart.f1 = 2 * bart.prec * bart.rec / (bart.prec + bart.rec)
      
      res$bart.precision = bart.prec
      res$bart.recall = bart.rec
      res$bart.f1 = bart.f1
    }
  }
  
  #------------------------------------------
  # plot
  if (plot) {
    categorical.idx = which(sapply(x.train, function(s) {is.factor(s)}))
    categorical.names = names(categorical.idx)
    
    if (num.var.plot == Inf | num.var.plot > ncol(x.train)){
      num.var.plot = ncol(x.train)
    }
    
    ## dart
    non.zero.idx = which(dart.pvip > 0)[1:min(num.var.plot, length(which(dart.pvip > 0)))]
    plot.n = length(non.zero.idx)
    if(length(non.zero.idx) < length(dart.pvip)) 
      warning(paste(length(which(dart.pvip == 0)), 
                    "predictors with marginal posterior variable inclusion probability of 0 omitted from plots."))
    
    plot(1:plot.n, dart.pvip[non.zero.idx], 
         type = "n", xlab = "Predictors", xaxt = "n", ylim = c(0, max(dart.pvip) * 1.1),
         ylab = "DART Posterior Variable Inclusion Probability")
    axis(1, at = 1:plot.n, labels = names(dart.pvip[non.zero.idx]), las = 2)
    for (j in non.zero.idx){
      points(j, dart.pvip[j],
             pch = ifelse(dart.pvip[j] <= 0.5, 1, 16),
             col = ifelse(names(dart.pvip[j]) %in% categorical.names, 'green', 'red'))
    }
    abline(h = 0.5, col = "grey", lty = 2)
    
    if (vip.selection) {
      ## bart
      non.zero.idx = which(bart.vip > 0)[1:min(num.var.plot, length(which(bart.vip > 0)))]
      plot.n = length(non.zero.idx)
      if(length(non.zero.idx) < length(bart.vip)) 
        warning(paste(length(which(bart.vip == 0)), 
                      "predictors with variable inclusion proprotion of 0 omitted from plots."))
      
      plot(1:plot.n, bart.vip[non.zero.idx], 
           type = "n", xlab = "Predictors", xaxt = "n", ylim = c(0, max(bart.vip) * 1.1),
           ylab = "BART VIP")
      axis(1, at = 1:plot.n, labels = names(bart.vip[non.zero.idx]), las = 2)
      for (j in non.zero.idx){
        points(j, bart.vip[j],
               pch = ifelse(bart.vip[j] <= (1/dim(x.train)[2]), 1, 16),
               col = ifelse(names(bart.vip[j]) %in% categorical.names, 'green', 'red'))
      }
      abline(h = 1 / dim(x.train)[2], col = "grey", lty = 2)
    }
  }
  
  #------------------------------
  # timer ends
  end = Sys.time()
  cat("Elapsed", end-start, '\n')
  
  return(res)
}