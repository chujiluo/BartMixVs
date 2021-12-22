mc.abc.vs = function(x, 
                     y, 
                     nabc=1000, 
                     tolerance=0.1, 
                     threshold=0.25,
                     beta.params=c(1.0, 1.0),
                     split.ratio=0.5, 
                     probit=F, 
                     true.idx=NULL,
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
                     verbose=FALSE,
                     mc.cores = 2L, 
                     nice = 19L, 
                     seed = 99L) {
  
  #------------------------------
  # timer starts
  start = Sys.time()

  res = list()
  
  #------------------------------
  # parallel setting
  if(.Platform$OS.type!='unix')
    stop('parallel::mcparallel/mccollect do not exist on windows')
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  parallel::mc.reset.stream()
  
  mc.cores.detected = detectCores()
  if(mc.cores > mc.cores.detected) mc.cores = mc.cores.detected
  if(nabc < mc.cores) mc.cores = nabc
  
  #------------------------------
  # parallel running
  beta.theta = rbeta(1, beta.params[1], beta.params[2])
  res$theta = beta.theta
  
  for(j in 1:mc.cores) {
    
    if (j <= nabc %% mc.cores){
      mc.nabc = floor(nabc / mc.cores) + 1
    }else{
      mc.nabc = floor(nabc / mc.cores)
    }
    
    parallel::mcparallel({psnice(value = nice);
      abc.vs(x = x, y = y, nabc = mc.nabc, tolerance = tolerance, threshold = threshold,
             beta.params = beta.params, beta.theta = beta.theta, split.ratio = split.ratio, 
             probit = probit, true.idx = true.idx, analysis = FALSE, sparse = sparse,
             xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, 
             rm.const = rm.const, k = k, power = power, base = base, split.prob = split.prob,
             ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery, 
             printevery = printevery, verbose = verbose)},
      silent = (j != 1))
  }
  
  results.list = parallel::mccollect()
  
  # collect parallel results
  models = results.list[[1]]$models
  actual.models = results.list[[1]]$actual.models
  model.errors = results.list[[1]]$model.errors
  
  if(mc.cores > 1) {
    for (j in 2:mc.cores) {
      models = rbind(models, results.list[[j]]$models)
      actual.models = rbind(actual.models, results.list[[j]]$actual.models)
      model.errors = c(model.errors, results.list[[j]]$model.errors)
    }
  }
  
  res$models = models
  res$actual.models = actual.models
  res$model.errors = model.errors
  
  rm(results.list)
  
  
  #------------------------------
  # top models
  ntop = ceiling(nabc * tolerance)
  idx0 = sort(model.errors, index.return = TRUE)$ix
  idx = idx0[1 : ntop]
  
  top.models = models[idx, ]
  top.actual.models = actual.models[idx, ]
  
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
  
  #------------------------------
  # timer
  end = Sys.time()
  cat("Elapsed", end-start, '\n')
  
  return(res)
}