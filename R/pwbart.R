pwbart = function(x.test,		         #x matrix to predict at
                  treedraws,		      #$treedraws from wbart or pbart
                  rm.const,            #$rm.const from wbart or pbart
                  mu=0,		            #mean to add on
                  mc.cores=1L,         #thread count
                  transposed=FALSE,	
                  dodraws=TRUE,
                  nice=19L             #mc.pwbart only
){
   if(!transposed) x.test = t(bartModelMatrix(x.test)[ , rm.const])
   
   p = length(treedraws$cutpoints)
   
   if(p != nrow(x.test))
      stop(paste0('The number of columns in x.test must be equal to ', p))
   
   res = .Call("cpwbart",
               treedraws,	#trees list
               x.test,      #the test x
               mc.cores   	#thread count
   )
   if(dodraws) return(res$yhat.test+mu)
   else return(apply(res$yhat.test, 2, mean)+mu)
}
