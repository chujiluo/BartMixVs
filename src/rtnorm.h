/*
 *  This file is copied from:
 *  
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  and modified by:
 *  Chuji Luo and Michael J. Daniels
 */
#ifndef GUARD_rtnorm
#define GUARD_rtnorm

#include "common.h"

double rtnorm(double mean, double tau, double sd, rn& gen);

RcppExport SEXP crtnorm(SEXP, SEXP, SEXP, SEXP);

#endif