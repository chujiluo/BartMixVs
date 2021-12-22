predict.wbart <- function(object, newdata, mc.cores=1, openmp=(mc.cores.openmp()>0), ...) {

    # p <- length(object$treedraws$cutpoints)
    # 
    # if(p!=ncol(newdata))
    #     stop(paste0('The number of columns in newdata must be equal to ', p))

    if(.Platform$OS.type == "unix") mc.cores.detected = detectCores()
    else mc.cores.detected = NA

    if(!is.na(mc.cores.detected) && mc.cores>mc.cores.detected) mc.cores = mc.cores.detected

    if(.Platform$OS.type != "unix" || openmp || mc.cores==1) call = pwbart
    else call = mc.pwbart

    if(length(object$mu) == 0) object$mu = object$offset

    return(call(newdata, object$treedraws, object$rm.const, mc.cores=mc.cores, mu=object$mu, ...))
}

