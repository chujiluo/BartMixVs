wbart = function(x.train, 
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
                 grp=NULL,    ## group indices for predictors, see bartModelMatrix() for details
                 xnames=NULL,  ## column names of x.train
                 categorical.idx=NULL,  ## column indices of categorical predictors in x.train
                 power=2.0,  ## not used when split.prob = "exponential"
                 base=-1.0,  ## by default, base = 0.95 when split.prob="polynomial" and base = 0.5 when split.prob = "exponential"
                 split.prob="polynomial",  ## specify if the split probability is polynomial (Chipman et al. 2010) or exponential (Rockova and Saha 2019)
                 k=2.0, 
                 sigmaf=NA,
                 sigest=NA, 
                 sigdf=3, 
                 sigquant=.90, 
                 lambda=NA,
                 fmean=mean(y.train), 
                 w=rep(1,length(y.train)),
                 ntree=200L, 
                 ndpost=1000L, 
                 nskip=1000L, 
                 keepevery=1L,
                 nkeeptrain=ndpost, 
                 nkeeptest=ndpost,
                 nkeeptestmean=ndpost, 
                 nkeeptreedraws=ndpost,
                 printevery=100L, 
                 transposed=FALSE,
                 verbose=TRUE  ## control if information is printed out
                 ) {
  #--------------------------------------------------
  #data
  n = length(y.train)
  
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
      p0 = length(unique(grp))  # number of predictors before dummification
    }
    categorical.idx = unique(grp)[which(sapply(unique(grp), function(s) sum(s==grp)) > 1)]
    
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
  if(length(rho) == 0) rho=p
  
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
  
  y.train = y.train - fmean
  
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
  if((nkeeptestmean != 0) & ((ndpost %% nkeeptestmean) != 0)) {
    nkeeptestmean = ndpost
    cat('*****nkeeptestmean set to ndpost\n')
  }
  if((nkeeptreedraws != 0) & ((ndpost %% nkeeptreedraws) != 0)) {
    nkeeptreedraws = ndpost
    cat('*****nkeeptreedraws set to ndpost\n')
  }
  
  #--------------------------------------------------
  #prior
  nu = sigdf
  if(is.na(lambda)) {
    if(is.na(sigest)) {
      if(p < n) {
        df = data.frame(t(x.train), y.train)
        lmf = lm(y.train~., df)
        sigest = summary(lmf)$sigma
      } else {
        sigest = sd(y.train)
      }
    }
    qchi = qchisq(1.0 - sigquant, nu)
    lambda = (sigest * sigest * qchi) / nu #lambda parameter for sigma prior
  } else {
    sigest = sqrt(lambda)
  }
  
  if(is.na(sigmaf)) {
    tau = (max(y.train) - min(y.train)) / (2 * k * sqrt(ntree))
  } else {
    tau = sigmaf / sqrt(ntree)
  }
  
  #--------------------------------------------------
  #call c++ function
  ptm = proc.time()
  res = .Call("cwbart",
              n,  #number of observations in training data
              p,  #dimension of x
              np, #number of observations in test data
              x.train,   #pxn training data x
              y.train,   #pxn training data x
              x.test,   #p*np test data x
              ntree,
              numcut,
              ndpost*keepevery,
              nskip,
              power,
              base,
              tau,
              nu,
              lambda,
              sigest,
              w,
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
              nkeeptestmean,
              nkeeptreedraws,
              printevery,
              xinfo,
              verbose
  )
  
  res$proc.time = proc.time() - ptm
  
  #--------------------------------------------------
  #returns
  res$mu = fmean
  res$yhat.train.mean = res$yhat.train.mean + fmean
  res$yhat.train = res$yhat.train + fmean
  res$yhat.test.mean = res$yhat.test.mean + fmean
  res$yhat.test = res$yhat.test + fmean
  
  if(nkeeptreedraws > 0)
    names(res$treedraws$cutpoints) = xnames
  
  #--------------------------------------------------
  #importance
  if(length(grp) == length(unique(grp))) {
    ## no dummy variables
    dimnames(res$varcount)[[2]] = as.list(xnames)
    dimnames(res$varprob)[[2]] = as.list(xnames)
    
    ## vip: variable inclusion proportions
    res$vip = colMeans(t(apply(res$varcount, 1, function(s) s / sum(s))))
    
    ## (marginal) posterior variable inclusion probability
    res$pvip = colMeans(res$varcount > 0)
    
    ## posterior s_j's (only in DART)
    res$varprob.mean = colMeans(res$varprob)
    
    ## mi: Metropolis importance
    mr.vecs = lapply(res$mr_vecs, function(s) lapply(s, function(v) v[-1]))  # remove the meaningless first 0
    res$mr_vecs = NULL
    res$mr.vecs = mr.vecs
    mr.mean = matrix(unlist(lapply(mr.vecs, function(s) lapply(s, function(v) ifelse(length(v) > 0, mean(v), 0.0)))), 
                     ncol = p, byrow = TRUE)
    res$mr.mean = mr.mean
    res$mi = colMeans(t(apply(mr.mean, 1, function(s) s / sum(s))))
    names(res$mi) = as.list(xnames)
    dimnames(res$mr.mean)[[2]] = as.list(xnames)
    
    ## untruncated mi
    # mr0_vecs = lapply(res$mr0_vecs, function(s) lapply(s, function(v) v[-1]))  # remove the meaningless first 0
    # res$mr0_vecs = mr0_vecs
    # mr0mean = matrix(unlist(lapply(mr0_vecs, function(s) lapply(s, function(v) ifelse(length(v)>0, mean(v), 0.0)))), 
    #                 ncol = p, byrow = TRUE)
    # res$mr0mean = mr0mean
    # res$mi0 = colMeans(t(apply(mr0mean, 1, function(s) s/sum(s))))
    # names(res$mi0) = as.list(xnames)
    # dimnames(res$mr0mean)[[2]] = as.list(xnames)
  } else {
    ## merge importance scores for dummy variables
    varcount = matrix(NA, nrow = nkeeptreedraws, ncol = p0)
    varprob = matrix(NA, nrow = nkeeptreedraws, ncol = p0)
    
    mr.vecs = lapply(res$mr_vecs, function(s) list(s[[1]][-1]))
    #mr0_vecs = lapply(res$mr0_vecs, function(s) list(s[[1]][-1]))
    varcount[, 1] = res$varcount[, 1]
    varprob[, 1] = res$varprob[, 1]
    
    j = 1
    for (l in 2:p) {
      if (grp[l] == grp[l-1]) {
        varcount[, j] = varcount[, j] + res$varcount[, l]
        varprob[, j] = varprob[, j] + res$varprob[, l]
        for (i in 1:nkeeptreedraws) {
          mr.vecs[[i]][[j]] = c(mr.vecs[[i]][[j]], res$mr_vecs[[i]][[l]][-1])
          #mr0_vecs[[i]][[j]] = c(mr0_vecs[[i]][[j]], res$mr0_vecs[[i]][[l]][-1])
        }
      } else {
        j = j + 1
        varcount[, j] = res$varcount[, l]
        varprob[, j] = res$varprob[, l]
        for (i in 1:nkeeptreedraws) {
          mr.vecs[[i]][[j]] = res$mr_vecs[[i]][[l]][-1]
          #mr0_vecs[[i]][[j]] = res$mr0_vecs[[i]][[l]][-1]
        }
      }
    }
    
    dimnames(varcount)[[2]] = as.list(xnames)
    dimnames(varprob)[[2]] = as.list(xnames)
    
    res$varcount = varcount
    res$varprob = varprob
    res$mr.vecs = mr.vecs
    res$mr_vecs = NULL
    #res$mr0_vecs = mr0_vecs
    
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
    
    ## untruncated mi 
    # mr0mean = matrix(unlist(lapply(mr0_vecs, function(s) lapply(s, function(v) ifelse(length(v)>0, mean(v), 0.0)))), 
    #                 ncol = p0, byrow = TRUE)
    # res$mr0mean = mr0mean
    # res$mi0 = colMeans(t(apply(mr0mean, 1, function(s) s/sum(s))))
    # dimnames(res$mr0mean)[[2]] = as.list(xnames)
    # names(res$mi0) = as.list(xnames)
  }
  
  res$rm.const = rm.const
  
  attr(res, 'class') = 'wbart'
  return(res)
}