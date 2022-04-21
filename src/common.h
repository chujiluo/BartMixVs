/*
 *  This file is copied from:
 *  
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  and modified by:
 *  Chuji Luo and Michael J. Daniels
 */
#ifndef GUARD_common_h
#define GUARD_common_h

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <ctime>
#include <cstddef>
#include <typeinfo>
#include <map>
#include <algorithm>

using std::endl;

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Rcpp.h>

#define printf Rprintf

#define cout Rcpp::Rcout

// log(2*pi)
#define LTPI 1.837877066409345483560659472811

// log(pi)
#define logpi 1.1447298858494001741434273513530587

// sqrt(2*pi)
#define RTPI 2.506628274631000502415765284811

#include "rn.h"

#endif