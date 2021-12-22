mc.pwbart = function(x.test,		#x matrix to predict at
                     treedraws,		#$treedraws from wbart
                     rm.const,      #$rm.const from wbart or pbart
                     mu=0,	     	#mean to add on
                     mc.cores=2L,
                     transposed=FALSE,
                     dodraws=TRUE,
                     nice=19L) {
    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')
    
    if(!transposed) x.test = t(bartModelMatrix(x.test)[ , rm.const])
    
    p = length(treedraws$cutpoints)
    
    if(p != nrow(x.test))
        stop(paste0('The number of columns in x.test must be equal to ', p))
    
    mc.cores.detected = detectCores()
    
    if(!is.na(mc.cores.detected) && (mc.cores > mc.cores.detected)) mc.cores = mc.cores.detected
    
    K = ncol(x.test)
    k = K%/%mc.cores - 1
    j = K
    for(i in 1:mc.cores) {
        if(i == mc.cores) h = 1
        else h = j - k
        
        parallel::mcparallel({psnice(value = nice);
            pwbart(matrix(x.test[ , h:j], nrow = p, ncol = j-h+1), treedraws, rm.const, mu, 1, TRUE)},
            silent = (i != 1))
        j = h - 1
    }
    
    pred.list = parallel::mccollect()
    
    pred = pred.list[[1]]
    
    if(mc.cores > 1) for(i in 2:mc.cores) pred = cbind(pred, pred.list[[i]])
    
    if(dodraws) return(pred+mu)
    else return(apply(pred, 2, mean) + mu)
}
