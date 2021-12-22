#ifndef GUARD_heterbart_h
#define GUARD_heterbart_h

#include "bart.h"
#include "heterbartfuns.h"
#include "heterbd.h"

class heterbart : public bart
{
  public:
    heterbart():bart() {}
    heterbart(size_t m):bart(m) {}
    void pr();
    void draw(double *sigma, rn& gen);
};

#endif