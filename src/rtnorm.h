#ifndef GUARD_rtnorm
#define GUARD_rtnorm

#include "common.h"

double rtnorm(double mean, double tau, double sd, rn& gen);

RcppExport SEXP crtnorm(SEXP, SEXP, SEXP, SEXP);

#endif