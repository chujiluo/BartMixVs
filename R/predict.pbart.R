predict.pbart <- function(object, newdata, mc.cores=1, openmp=(mc.cores.openmp()>0), ...) {

    # p <- length(object$treedraws$cutpoints)
    # 
    # if(p!=ncol(newdata))
    #     stop(paste0('The number of columns in newdata must be equal to ', p))

    if(.Platform$OS.type == "unix") mc.cores.detected = detectCores()
    else mc.cores.detected = NA

    if(!is.na(mc.cores.detected) && mc.cores>mc.cores.detected) mc.cores = mc.cores.detected

    if(.Platform$OS.type != "unix" || openmp || mc.cores==1) call = pwbart
    else call = mc.pwbart

    if(length(object$binaryOffset)==0) object$binaryOffset=object$offset

    pred = list(yhat.test=call(newdata, object$treedraws, object$rm.const, mc.cores=mc.cores,
                                mu=object$binaryOffset, ...))

    pred$prob.test = pnorm(pred$yhat.test)
    pred$prob.test.mean = apply(pred$prob.test, 2, mean)
    pred$binaryOffset = object$binaryOffset
    attr(pred, 'class') = 'pbart'

    return(pred)
}

