abc.vs = function(x, 
                  y, 
                  nabc=1000, 
                  tolerance=0.1, 
                  threshold=0.25,
                  beta.params=c(1.0, 1.0), 
                  beta.theta=NA,
                  split.ratio=0.5, 
                  probit=F, 
                  true.idx=NULL,
                  analysis=TRUE,
                  sparse = FALSE,
                  xinfo=matrix(0.0,0,0), 
                  numcut=100L,
                  usequants=FALSE, 
                  cont=FALSE, 
                  rm.const=TRUE, 
                  k=2.0, 
                  power=2.0, 
                  base=0.95,
                  split.prob="polynomial",
                  ntree=10L, 
                  ndpost=1, 
                  nskip=200,
                  keepevery=1L, 
                  printevery=100L, 
                  verbose=FALSE) {
  
  res = list()
  
  #------------------------------
  # timer starts
  #start = Sys.time()
  
  #------------------------------
  # data
  n = nrow(x)
  p = ncol(x)
  
  #------------------------------
  # returns
  models = matrix(NA, nrow = 1, ncol = p)  ## pool of available predictors
  actual.models = matrix(NA, nrow = 1, ncol = p)  ## predictors used in the BART ensemble
  model.errors = c()  ## MSE for continuous y and mean log loss for binary y
  
  #------------------------------
  # theta ~ Beta(betaparams)
  if(is.na(beta.theta)) beta.theta = rbeta(1, beta.params[1], beta.params[2])
  res$theta = beta.theta
  
  #------------------------------
  # run ABC Bayesian forest
  for (i in 1:nabc) {
    ##---------------------
    ## split data
    train.idx = sample(1:n, size = floor(split.ratio * n), replace = FALSE)
    x.train = x[train.idx, ]
    y.train = y[train.idx]
    n.train = length(y.train)
    x.test = x[-train.idx, ]
    y.test = y[-train.idx]
    n.test = length(y.test)
    
    ##---------------------
    ## pick a subset
    current.model = c()
    
    ## make sure the model selected contains at least one predictor
    while (!any(current.model == 1)) {
      current.model = c()
      
      for (j in 1:p) {
        if(rbinom(1, 1, beta.theta))
          current.model[j] = 1
        else
          current.model[j] = 0
      }
    }
    
    if(i == 1)
      models[1, ] = current.model
    else
      models = rbind(models, current.model)
    
    ##---------------------
    ## bart
    if(probit) {
      bart.obj = pbart(x.train = x.train[, which(current.model == 1)], y.train = y.train, 
                       x.test = x.test[, which(current.model == 1)], sparse = sparse,
                       xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, rm.const = rm.const, 
                       k = k, power = power, base = base, split.prob = split.prob,
                       ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, 
                       printevery = printevery, verbose = verbose)
      
      ## model error: mean logarithmic loss
      prob.test = bart.obj$prob.test
      model.error = 0
      for (j in 1:n.test) {
        if(((prob.test[j] == 0) & (y.test[j] == 0)) | ((prob.test[j] == 1) & (y.test[j] == 1))){
          model.error = model.error
        } else {
          model.error = model.error + y.test[j] * log(prob.test[j]) + (1 - y.test[j]) * log(1 - prob.test[j])
        }
      }
      model.errors[i] = - model.error / n.test
      
      actual.current.model = rep(0, p)
      actual.current.model[which(current.model == 1)[which(bart.obj$varcount > 0)]] = 1
      
      if(i == 1)
        actual.models[1, ] = actual.current.model
      else
        actual.models = rbind(actual.models, actual.current.model)
      
    }else{
      bart.obj = wbart(x.train = x.train[, which(current.model == 1)], y.train = y.train, 
                       x.test = x.test[, which(current.model == 1)], sparse = sparse,
                       xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, rm.const = rm.const, 
                       k = k, power = power, base = base, split.prob = split.prob,
                       ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, 
                       printevery = printevery, verbose = verbose)
      
      pseudo.y.test = bart.obj$yhat.test + rnorm(n.test, 0, bart.obj$sigma[nskip + ndpost])
      model.errors[i] = sqrt(mean((pseudo.y.test - y.test) ^ 2))
      
      actual.current.model = rep(0, p)
      actual.current.model[which(current.model==1)[which(bart.obj$varcount > 0)]] = 1
      
      if(i == 1)
        actual.models[1, ] = actual.current.model
      else
        actual.models = rbind(actual.models, actual.current.model)
    }
    
  }
  
  if(analysis) {
    #------------------------------
    # top models
    ntop = ceiling(nabc * tolerance)
    idx0 = sort(model.errors, index.return = TRUE)$ix
    idx = idx0[1 : ntop]
    
    top.models = models[idx, ]
    top.actual.models = actual.models[idx, ]
    
    res$models = models
    res$actual.models = actual.models
    res$model.errors = model.errors
    res$idx = idx
    res$top.models = top.models
    res$top.actual.models = top.actual.models
    
    #------------------------------
    # marginal inclusion probabilities
    mip = colMeans(top.actual.models)
    best.model = which(mip >= threshold)
    
    res$mip = mip
    res$best.model = best.model
    
    #------------------------------
    # score
    if(length(true.idx) > 0) {
      true.len = length(true.idx)
      
      tp = length(which(best.model %in% true.idx))
      positive.len = length(best.model)
      
      res$precision = (tp * 1.0) / (positive.len * 1.0)
      res$recall = (tp * 1.0) / (true.len * 1.0)
      res$f1 = 2 * res$precision * res$recall / (res$precision + res$recall)
    }
  }else{
    res$models = models
    res$actual.models = actual.models
    res$model.errors = model.errors
  }
 
  
  #------------------------------
  # timer ends
  #end = Sys.time()
  #cat("Elapsed", end-start, '\n')
  
  return(res)
}