#ifndef GUARD_heterbd_h
#define GUARD_heterbd_h

#include "info.h"
#include "tree.h"
#include "treefuns.h"
#include "bartfuns.h"
#include "heterbartfuns.h"

bool heterbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, 
             std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen,
             std::vector<std::vector<double>>& mr_vec);

#endif