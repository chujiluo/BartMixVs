/*
 *  This file is copied from:
 *  
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  and modified by:
 *  Chuji Luo and Michael J. Daniels
 */
#ifndef GUARD_bd_h
#define GUARD_bd_h

#include "info.h"
#include "tree.h"
#include "treefuns.h"
#include "bartfuns.h"

bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma, 
        std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen, 
        std::vector<std::vector<double>>& mr_vec);

#endif