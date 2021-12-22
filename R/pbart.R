pbart = function(x.train, 
                 y.train, 
                 x.test=matrix(0.0,0,0),
                 sparse=FALSE, 
                 theta=0, 
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
                 grp = NULL, 
                 xnames = NULL, 
                 categorical.idx = NULL,
                 k=2.0, 
                 power=2.0, 
                 base=-1.0, 
                 split.prob = "polynomial",
                 binaryOffset=NULL,
                 ntree=50L, 
                 ndpost=1000L, 
                 nskip=1000L, 
                 keepevery=1L,
                 nkeeptrain=ndpost, 
                 nkeeptest=ndpost, 
                 nkeeptreedraws=ndpost,
                 printevery=100L, 
                 transposed=FALSE,
                 verbose=TRUE) {
  #--------------------------------------------------
  #data
  n = length(y.train)
  
  if(length(binaryOffset) == 0) binaryOffset = qnorm(mean(y.train))
  
  if(!transposed) {
    if(is.null(dim(x.train))) {
      xnames = "X"
    } else {
      xnames = dimnames(x.train)[[2]]
    }
    
    temp = bartModelMatrix(x.train, numcut, usequants = usequants,
                           cont = cont, xinfo = xinfo, rm.const = rm.const)
    x.train = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    if(length(x.test) > 0) {
      x.test = bartModelMatrix(x.test)
      x.test = t(x.test[ , temp$rm.const])
    }
    rm.const = temp$rm.const
    grp = temp$grp
    rm(temp)
    
    if(length(grp) == 0){
      p0 = nrow(x.train)
      grp = 1:p0
    } else {
      p0 = length(unique(grp))
    }
    categorical.idx = unique(grp)[which(sapply(unique(grp), function(s) sum(s == grp)) > 1)]
  
  } else {
    if(any(length(rm.const) == 0, length(grp) == 0, length(xnames) == 0))
      stop('Did not provide rm.const, grp and xnames for x.train after transpose!')
    if(is.logical(rm.const))
      stop('Did not provide rm.const for x.train after transpose!')
    if((length(grp) > length(unique(grp))) & (length(categorical.idx) <= 0))
      stop('Did not provide categorical.idx for x.train that contains categorical predictors!')
    
    p0 = length(unique(grp))
  }
  
  if(n != ncol(x.train))
    stop('The length of y.train and the number of rows in x.train must be identical')
  
  p = nrow(x.train)
  np = ncol(x.test)
  if(length(rho) == 0) rho <- p
  
  if(!(split.prob %in% c("polynomial", "exponential"))) {
    stop("split.prob is either polynomial or exponential.")
  } else {
    if(split.prob == "polynomial") {
      if(base < 0)
        base = 0.95
    }
    if(split.prob == "exponential") {
      power = -1.0
      if(base < 0)
        base = 0.5
    }
  }
  
  #--------------------------------------------------
  #set nkeeps for thinning
  if((nkeeptrain != 0) & ((ndpost %% nkeeptrain) != 0)) {
    nkeeptrain = ndpost
    cat('*****nkeeptrain set to ndpost\n')
  }
  if((nkeeptest != 0) & ((ndpost %% nkeeptest) != 0)) {
    nkeeptest = ndpost
    cat('*****nkeeptest set to ndpost\n')
  }
  if((nkeeptreedraws != 0) & ((ndpost %% nkeeptreedraws) != 0)) {
    nkeeptreedraws = ndpost
    cat('*****nkeeptreedraws set to ndpost\n')
  }
  
  
  #--------------------------------------------------
  ptm = proc.time()
  #call
  res = .Call("cpbart",
              n,  #number of observations in training data
              p,  #dimension of x
              np, #number of observations in test data
              x.train,   #p*n training data x
              y.train,   #n*1 training data y
              x.test,    #p*np test data x
              ntree,
              numcut,
              ndpost*keepevery,
              nskip,
              power,
              base,
              binaryOffset,
              3/(k*sqrt(ntree)),
              sparse,
              theta,
              omega,
              grp,
              a,
              b,
              rho,
              augment,
              nkeeptrain,
              nkeeptest,
              nkeeptreedraws,
              printevery,
              xinfo,
              verbose
  )
  
  res$proc.time = proc.time() - ptm
  
  #---------------------------------------------------------
  #returns
  if(nkeeptrain > 0) {
    res$yhat.train = res$yhat.train + binaryOffset
    res$prob.train = pnorm(res$yhat.train)
    res$prob.train.mean = apply(res$prob.train, 2, mean)
  } else {
    res$yhat.train = NULL
  }
  
  if(np > 0) {
    res$yhat.test = res$yhat.test + binaryOffset
    res$prob.test = pnorm(res$yhat.test)
    res$prob.test.mean = apply(res$prob.test, 2, mean)
  } else {
    res$yhat.test = NULL
  }
  
  if(nkeeptreedraws > 0)
    names(res$treedraws$cutpoints) = xnames
  
  #--------------------------------------------------
  #importance
  if(length(grp) == length(unique(grp))) {
    ## no dummy variables
    dimnames(res$varcount)[[2]] = as.list(xnames)
    dimnames(res$varprob)[[2]] = as.list(xnames)
    
    ## vip
    res$vip = colMeans(t(apply(res$varcount, 1, function(s) s / sum(s))))
    
    ## (marginal) posterior variable inclusion probability
    res$pvip = colMeans(res$varcount > 0)
    
    ## posterior s_j's (in DART)
    res$varprob.mean = colMeans(res$varprob)
    
    ## mi
    mr.vecs = lapply(res$mr_vecs, function(s) lapply(s, function(v) v[-1]))  # remove the meaningless first 0
    res$mr.vecs = mr.vecs
    res$mr_vecs = NULL
    mr.mean = matrix(unlist(lapply(mr.vecs, function(s) lapply(s, function(v) ifelse(length(v) > 0, mean(v), 0.0)))), 
                     ncol = p, byrow = TRUE)
    res$mr.mean = mr.mean
    res$mi = colMeans(t(apply(mr.mean, 1, function(s) s / sum(s))))
    names(res$mi) = as.list(xnames)
    dimnames(res$mr.mean)[[2]] = as.list(xnames)
  } else {
    ## merge importance scores for dummy variables
    varcount = matrix(NA, nrow = nkeeptreedraws, ncol = p0)
    varprob = matrix(NA, nrow = nkeeptreedraws, ncol = p0)
    
    mr.vecs = lapply(res$mr_vecs, function(s) list(s[[1]][-1]))
    varcount[, 1] = res$varcount[, 1]
    varprob[, 1] = res$varprob[, 1]
    
    j = 1
    for (l in 2:p) {
      if (grp[l] == grp[l-1]) {
        varcount[, j] = varcount[, j] + res$varcount[, l]
        varprob[, j] = varprob[, j] + res$varprob[, l]
        for (i in 1:nkeeptreedraws) {
          mr.vecs[[i]][[j]] = c(mr.vecs[[i]][[j]], res$mr_vecs[[i]][[l]][-1])
        }
      } else {
        j = j + 1
        varcount[, j] = res$varcount[, l]
        varprob[, j] = res$varprob[, l]
        for (i in 1:nkeeptreedraws) {
          mr.vecs[[i]][[j]] = res$mr_vecs[[i]][[l]][-1]
        }
      }
    }
    
    dimnames(varcount)[[2]] = as.list(xnames)
    dimnames(varprob)[[2]] = as.list(xnames)
    
    res$varcount = varcount
    res$varprob = varprob
    res$mr.vecs = mr.vecs
    res$mr_vecs = NULL
    
    ## vip
    res$vip = colMeans(t(apply(varcount, 1, function(s) s / sum(s))))
    
    ## within-type vip
    within.type.vip = rep(0, p0)
    for (i in 1:nkeeptreedraws) {
      if (sum(varcount[i, categorical.idx]) != 0) {
        within.type.vip[categorical.idx] = within.type.vip[categorical.idx] + 
          varcount[i, categorical.idx] / sum(varcount[i, categorical.idx])
      }
      if (sum(varcount[i, -categorical.idx]) != 0) {
        within.type.vip[-categorical.idx] = within.type.vip[-categorical.idx] + 
          varcount[i, -categorical.idx] / sum(varcount[i, -categorical.idx])
      }
    }
    res$within.type.vip = within.type.vip / nkeeptreedraws
    names(res$within.type.vip) = xnames
    
    ## (marginal) posterior variable inclusion probability
    res$pvip = colMeans(varcount > 0)
    
    ## posterior s_j's (in DART)
    res$varprob.mean = colMeans(varprob)
    
    ## mi
    mr.mean = matrix(unlist(lapply(mr.vecs, function(s) lapply(s, function(v) ifelse(length(v) > 0, mean(v), 0.0)))), 
                     ncol = p0, byrow = TRUE)
    res$mr.mean = mr.mean
    res$mi = colMeans(t(apply(mr.mean, 1, function(s) s / sum(s))))
    dimnames(res$mr.mean)[[2]] = as.list(xnames)
    names(res$mi) = as.list(xnames)
  }
  
  
  res$rm.const = rm.const
  res$binaryOffset = binaryOffset
  
  attr(res, 'class') = 'pbart'
  return(res)
}