mc.pbart = function(x.train, 
                    y.train, 
                    x.test = matrix(0.0, 0L, 0L),
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
                    k=2.0, ## BEWARE: do NOT use k for other purposes below
                    power=2.0, 
                    base=0.95,
                    split.prob="polynomial",
                    binaryOffset=NULL,
                    ntree=50L, 
                    ndpost=1000L, 
                    nskip=100L,
                    keepevery=1L, 
                    printevery=100L,
                    keeptrainfits=TRUE, 
                    transposed=FALSE,
                    verbose=FALSE,
                    mc.cores = 2L, 
                    nice = 19L,
                    seed = 99L) {
    
    # data
    if(is.null(dim(x.train))) {
        xnames = "X"
    } else {
        xnames = dimnames(x.train)[[2]]  # predictor names before dummification
    }
    
    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')
    
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()
    
    if(!transposed) {
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
        
        p0 = length(unique(grp))  # number of predictors before dummification
        categorical.idx = unique(grp)[which(sapply(unique(grp), 
                                                   function(s) sum(s == grp)) > 1)]
        
    }
    
    mc.cores.detected = detectCores()
    
    if(mc.cores > mc.cores.detected) mc.cores = mc.cores.detected
    
    mc.ndpost = ceiling(ndpost/mc.cores)
    
    for(i in 1:mc.cores) {
        parallel::mcparallel({psnice(value = nice);
            pbart(x.train = x.train, y.train = y.train, x.test = x.test,
                  sparse = sparse, theta = theta, omega = omega, a = a, b = b, augment = augment, rho = rho,
                  xinfo = xinfo, numcut = numcut, rm.const = rm.const, 
                  grp = grp, xnames = xnames, categorical.idx = categorical.idx,
                  k = k, power = power, base = base, split.prob = split.prob, binaryOffset = binaryOffset,
                  ntree = ntree, ndpost = mc.ndpost, nskip = nskip, keepevery = keepevery,
                  printevery = printevery, transposed = TRUE, verbose = verbose)},
            silent = (i != 1))
        ## to avoid duplication of output
        ## capture stdout from first posterior only
    }
    
    post.list = parallel::mccollect()
    
    post = post.list[[1]]
    
    if(mc.cores == 1 | attr(post, 'class') != 'pbart') return(post)
    else {
        if(class(rm.const)[1]!='logical') post$rm.const = rm.const
        
        post$ndpost = mc.cores * mc.ndpost
        
        p = nrow(x.train[post$rm.const, ])
        
        old.text = paste0(as.character(mc.ndpost), ' ', as.character(ntree), ' ', as.character(p))
        old.stop = nchar(old.text)
        
        post$treedraws$trees = sub(old.text,
                                   paste0(as.character(post$ndpost), ' ',
                                          as.character(ntree), ' ',
                                          as.character(p)),
                                   post$treedraws$trees)
        
        post$vip = (post$vip) * mc.ndpost
        post$pvip = (post$pvip) * mc.ndpost
        post$varprob.mean = (post$varprob.mean) * mc.ndpost
        post$mi = (post$mi) * mc.ndpost
        if (length(grp) > length(unique(grp)))
            post$within.type.vip = (post$within.type.vip) * mc.ndpost
        
        keeptestfits = length(x.test)>0
        
        for(i in 2:mc.cores) {
            if(keeptrainfits) {
                post$yhat.train = rbind(post$yhat.train, post.list[[i]]$yhat.train)
                post$prob.train = rbind(post$prob.train, post.list[[i]]$prob.train)
            }
            
            if(keeptestfits) {
                post$yhat.test = rbind(post$yhat.test, post.list[[i]]$yhat.test)
                post$prob.test = rbind(post$prob.test, post.list[[i]]$prob.test)
            }
            
            post$varcount = rbind(post$varcount, post.list[[i]]$varcount)
            post$varprob = rbind(post$varprob, post.list[[i]]$varprob)
            
            post$mr.vecs = c(post$mr.vecs, post.list[[i]]$mr.vecs)
            post$mr.mean = rbind(post$mr.mean, post.list[[i]]$mr.mean)
            
            post$vip = post$vip + (post.list[[i]]$vip) * mc.ndpost
            post$pvip = post$pvip + (post.list[[i]]$pvip) * mc.ndpost
            post$varprob.mean = post$varprob.mean + (post.list[[i]]$varprob.mean) * mc.ndpost
            post$mi = post$mi + (post.list[[i]]$mi) * mc.ndpost
            if (length(grp) > length(unique(grp)))
                post$within.type.vip = post$within.type.vip + (post.list[[i]]$within.type.vip) * mc.ndpost
            
            post$treedraws$trees = paste0(post$treedraws$trees,
                                          substr(post.list[[i]]$treedraws$trees, old.stop + 2,
                                                 nchar(post.list[[i]]$treedraws$trees)))
        }
        
        if(length(post$prob.train.mean) > 0)
            post$prob.train.mean = apply(post$prob.train, 2, mean)
        
        if(length(post$prob.test.mean)>0)
            post$prob.test.mean = apply(post$prob.test, 2, mean)
        
        # process importance
        dimnames(post$varcount)[[2]] = as.list(xnames)
        dimnames(post$varprob)[[2]] = as.list(xnames)
        dimnames(post$mr.mean)[[2]] = as.list(xnames)
        
        post$vip = post$vip / (mc.ndpost * mc.cores)
        post$pvip = post$pvip / (mc.ndpost * mc.cores)
        post$varprob.mean = post$varprob.mean / (mc.ndpost * mc.cores)
        post$mi = post$mi / (mc.ndpost * mc.cores)
        if (length(grp) > length(unique(grp)))
            post$within.type.vip = post$within.type.vip / (mc.ndpost * mc.cores)
        
        names(post$vip) = as.list(xnames)
        names(post$pvip) = as.list(xnames)
        names(post$varprob.mean) = as.list(xnames)
        names(post$mi) = as.list(xnames)
        if (length(grp) > length(unique(grp)))
            names(post$within.type.vip) = as.list(xnames)
        
        attr(post, 'class') = 'pbart'
        
        return(post)
    }
}
