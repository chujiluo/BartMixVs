single.model.evaluation = function(x.train, 
                                   y.train, 
                                   x.test, 
                                   y.test, 
                                   n.train, 
                                   n.test, 
                                   n, 
                                   probit,
                                   xinfo, 
                                   numcut, 
                                   usequants, 
                                   cont, 
                                   rm.const, 
                                   k, 
                                   power, 
                                   base,
                                   split.prob,
                                   ntree, 
                                   ndpost, 
                                   nskip, 
                                   keepevery, 
                                   printevery,
                                   verbose) {
  if(probit) {
    bart = pbart(x.train = x.train, y.train = y.train, x.test = x.test, sparse = FALSE,
                 xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, rm.const = rm.const, 
                 k = k, power = power, base = base, split.prob = split.prob,
                 ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, 
                 printevery = printevery, verbose = verbose)
    
    ## model error: mean logarithmic loss
    prob.test = bart$prob.test.mean
    model.error = 0
    for (i in 1:n.test) {
      if((prob.test[i] == 0 & y.test[i] == 0) | (prob.test[i] == 1 & y.test[i] == 1)){
        model.error = model.error
      } else {
        model.error = model.error + y.test[i] * log(prob.test[i]) + (1 - y.test[i]) * log(1 - prob.test[i])
      }
    }
    model.error = - model.error / n.test
  } else {
    bart = wbart(x.train = x.train, y.train = y.train, x.test = x.test, sparse = FALSE,
                 xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, rm.const = rm.const, 
                 k = k, power = power, base = base, split.prob = split.prob,
                 ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, 
                 printevery = printevery, verbose = verbose)
    
    ## model error
    model.error = mean((bart$yhat.test.mean - y.test)^2)
  }
  
  # elpd.loo0
  loo.result = bart.loo(object = bart, y.train = y.train, probit = probit, n.train = n.train, 
                        nskip = nskip, ndpost = ndpost, keepevery = keepevery)
  elpd.loo = loo.result$estimates[1, 1]
  badK = sum(loo.result$diagnostics$pareto_k >= 0.7) / n.train
  
  res = list()
  res$model.error = model.error
  res$elpd.loo = elpd.loo
  res$badK = badK
  
  return(res)
}

multiple.models.evaluation = function(x.train, 
                                      y.train, 
                                      x.test, 
                                      y.test, 
                                      models,
                                      n.train, 
                                      n.test, 
                                      n, 
                                      probit,
                                      xinfo, 
                                      numcut, 
                                      usequants, 
                                      cont, 
                                      rm.const, 
                                      k, 
                                      power, 
                                      base,
                                      split.prob,
                                      ntree, 
                                      ndpost, 
                                      nskip, 
                                      keepevery, 
                                      printevery,
                                      verbose) {
  
  if(is.vector(models)) {
    res = matrix(NA, nrow = 1, ncol = 3)
    models = matrix(models, nrow = 1, ncol = length(models))
  } else {
    res = matrix(NA, nrow = nrow(models), ncol = 3)
  }
    
  for (i in 1:nrow(models)) {
    if(probit) {
      bart = pbart(x.train = x.train[, models[i, ], drop = FALSE], y.train = y.train, 
                   x.test = x.test[, models[i, ], drop = FALSE], 
                   sparse = FALSE, xinfo = xinfo, numcut = numcut, usequants = usequants, cont =cont, 
                   rm.const = rm.const, k = k, power = power, base = base, split.prob = split.prob,
                   ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, 
                   printevery = printevery, verbose = verbose)
      
      ## model error: mean logarithmic loss
      prob.test = bart$prob.test.mean
      res[i, 1] = 0
      for (j in 1:n.test) {
        if((prob.test[j] == 0 & y.test[j] == 0) | (prob.test[j] == 1 & y.test[j] == 1)){
          res[i, 1] = res[i, 1]
        } else {
          res[i, 1] = res[i, 1] + y.test[j] * log(prob.test[j]) + (1 - y.test[j]) * log(1 - prob.test[j])
        }
      }
      res[i, 1] = - res[i, 1] / n.test
    } else {
      bart = wbart(x.train = x.train[, models[i, ], drop = FALSE], y.train = y.train, 
                   x.test = x.test[, models[i, ], drop = FALSE], 
                   sparse = FALSE, xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, 
                   rm.const = rm.const, k = k, power = power, base = base, split.prob = split.prob,
                   ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, 
                   printevery = printevery, verbose = verbose)
      
      ## model error
      res[i, 1] = mean((bart$yhat.test.mean - y.test) ^ 2)
    }
    
    # elpd.loo0
    loo.result = bart.loo(object = bart, y.train = y.train, probit = probit, n.train = n.train, 
                          nskip = nskip, ndpost = ndpost, keepevery = keepevery)
    res[i, 2] = loo.result$estimates[1, 1]
    res[i, 3] = sum(loo.result$diagnostics$pareto_k >= 0.7) / n.train
  }
  
  result = list()
  result$res = res
  result$models = models
  return(result)
}

mc.backward.vs = function(x, 
                          y, 
                          split.ratio=0.8,
                          probit=F, 
                          true.idx=NULL,
                          xinfo=matrix(0.0,0,0), 
                          numcut=100L,
                          usequants=FALSE, 
                          cont=FALSE, 
                          rm.const=TRUE, 
                          k=2.0, 
                          power=2.0, 
                          base=0.95,
                          split.prob="polynomial",
                          ntree=50L, 
                          ndpost=1000, 
                          nskip=1000,
                          keepevery=1L,
                          printevery=100L, 
                          verbose=FALSE,
                          mc.cores = 2L, 
                          nice = 19L, 
                          seed = 99L) {
  
  #------------------------------
  # timer starts
  start = Sys.time()
  
  #------------------------------
  # parallel setting
  if(.Platform$OS.type!='unix')
    stop('parallel::mcparallel/mccollect do not exist on windows')
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  parallel::mc.reset.stream()
  
  mc.cores.detected = detectCores()
  if(mc.cores > mc.cores.detected) mc.cores = mc.cores.detected
  
  #------------------------------
  # data
  n = nrow(x)
  p = ncol(x)

  ## split data into training and test
  train.idx = sample(1:n, size = floor(split.ratio * n), replace = FALSE)
  x.train = x[train.idx, ]
  y.train = y[train.idx]
  n.train = length(y.train)
  x.test = x[-train.idx, ]
  y.test = y[-train.idx]
  n.test = length(y.test)

  #------------------------------
  # returns
  models = list()        # winner models from each iteration of backward, length = p
  model.errors = c()     # winner models' model errors, length = p
  elpd.loos = c()        # Bayesian loo of winner models from each iteration of backward, length = p
  badKs = c()            # badKs from Bayesian loo for winner models from each iteration of backward, length = p
  
  all.models = list()    # all evaluated models
  all.model.errors = c() # model errors corresponding to all.models (same order)

  #------------------------------
  # model 0: full model
  models[[1]] = 1:p
  all.models[[1]] = 1:p
  temp.result = single.model.evaluation(x.train = x.train, y.train = y.train, x.test = x.test, y.test = y.test, 
                                        n.train = n.train, n.test = n.test, n = n, probit = probit,
                                        xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, 
                                        rm.const = rm.const, 
                                        k = k, power = power, base = base, split.prob = split.prob,
                                        ntree = ntree, ndpost = ndpost, nskip = nskip, 
                                        keepevery = keepevery, printevery = printevery, verbose = verbose)
  model.errors[1] = temp.result$model.error
  all.model.errors[1] = temp.result$model.error
  elpd.loos[1] = temp.result$elpd.loo
  badKs[1] = temp.result$badK
  
  rm(temp.result)

  #------------------------------
  # backward in sequential, parallel within each iteration
  for (i in 1:(p-1)) {
    ## model[[i]]: winner model from last iteration
    cat("iteration =", (i+1), "/")
    
    num.models = length(models[[i]])
    reduced.models = matrix(NA, nrow = num.models, ncol = (num.models - 1))   # each row is a new model
    for (j in 1:num.models) {
      reduced.models[j, ] = models[[i]][-j]
    }
    
    ## parallel
    if (num.models < mc.cores) {
      mc.cores = num.models
      mc.num.models = 1
    } else {
      mc.num.models = floor(num.models / mc.cores)
    }
    
    for(j in 1:mc.cores) {
      if (j <= num.models %% mc.cores){
        start.row = (j - 1) * (mc.num.models + 1) + 1
        end.row = j * (mc.num.models + 1)
      }
      else{
        start.row = (j - 1) * mc.num.models + 1 + num.models %% mc.cores
        end.row = j * mc.num.models + num.models %% mc.cores
      }
      
      
      parallel::mcparallel({psnice(value = nice);
        multiple.models.evaluation(x.train = x.train, y.train = y.train, x.test = x.test, y.test = y.test, 
                                   models = reduced.models[start.row:end.row, , drop=FALSE],
                                   n.train = n.train, n.test = n.test, n = n, probit = probit,
                                   xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, 
                                   rm.const = rm.const, 
                                   k = k, power = power, base = base, split.prob = split.prob,
                                   ntree = ntree, ndpost = ndpost, nskip = nskip, 
                                   keepevery = keepevery, printevery = printevery, verbose = verbose)},
        silent = (j != 1))
    }
    
    reduced.results.list = parallel::mccollect()
    
    ## collect parallel results
    reduced.results = reduced.results.list[[1]]$res
    reduced.models = reduced.results.list[[1]]$models
    
    if(mc.cores > 1) {
      for (j in 2:mc.cores) {
        reduced.results = rbind(reduced.results, reduced.results.list[[j]]$res)
        reduced.models = rbind(reduced.models, reduced.results.list[[j]]$models)
      }
    }
    
    rm(reduced.results.list)
    
    ## pick the winner model at this iteration
    l.star = which.min(reduced.results[, 1])
    models[[i+1]] = reduced.models[l.star, ]
    model.errors[i+1] = reduced.results[l.star, 1]
    elpd.loos[i+1] = reduced.results[l.star, 2]
    badKs[i+1] = reduced.results[l.star, 3]

    ## store results of all evaluated models
    for (j in 1:num.models) {
      all.models = c(all.models, list(reduced.models[j, ]))
    }
    all.model.errors = c(all.model.errors, reduced.results[, 1])
    
    rm(reduced.results)
  }
  
  cat("\n")
  
  #------------------------------
  # pick the best model
  best.model.order = which.max(elpd.loos)
  best.model.cols = models[[best.model.order]]
  best.model.names = dimnames(x)[[2]][best.model.cols]
  
  #------------------------------
  # returns
  res = list()
  res$best.model.names = best.model.names
  res$best.model.cols = best.model.cols
  res$best.model.order = best.model.order
  res$models = models
  res$model.errors = model.errors
  res$elpd.loos = elpd.loos
  res$badKs = badKs
  res$all.models = all.models
  res$all.model.errors = all.model.errors

  #------------------------------
  # score the results
  if(length(true.idx) > 0) {
    true.len = length(true.idx)
    
    tp = length(which(best.model.cols %in% true.idx))
    positive.len = length(best.model.cols)
    
    res$precision = (tp * 1.0) / (positive.len * 1.0)
    res$recall = (tp * 1.0) / (true.len * 1.0)
    res$f1 = 2 * res$precision * res$recall / (res$precision + res$recall)
  }
  
  #------------------------------
  # distinguish acceptable and unacceptable models
  if(length(true.idx) > 0) {
    all.models.idx = c()    # true if acceptable; false if unacceptable
    for (i in 1:length(all.models)) {
      all.models.idx[i] = all(true.idx %in% all.models[[i]])
    }
    res$all.models.idx = all.models.idx
  }
  
  #------------------------------
  # timer
  end = Sys.time()
  cat("Elapsed", end-start, '\n')
  
  return(res)
}

