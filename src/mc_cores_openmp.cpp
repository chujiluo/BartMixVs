/*
 *  This file is copied from:
 *  
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  and modified by:
 *  Chuji Luo and Michael J. Daniels
 */
#include "common.h"

RcppExport SEXP mc_cores_openmp() {

#ifdef _OPENMP
    
    int mc_cores_openmp=omp_get_num_threads();
    
#else
    
    int mc_cores_openmp=0;
    
#endif
    
    return Rcpp::wrap(mc_cores_openmp);
    
}