#include "common.h"

RcppExport SEXP mc_cores_openmp() {

#ifdef _OPENMP
    
    int mc_cores_openmp=omp_get_num_threads();
    
#else
    
    int mc_cores_openmp=0;
    
#endif
    
    return Rcpp::wrap(mc_cores_openmp);
    
}