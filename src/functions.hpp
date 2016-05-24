#define ARMA_USE_LAPACK
#define ARMA_NO_DEBUG

// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <bigmemory/BigMatrix.h>

using namespace Rcpp;
using namespace arma;

#ifndef __FUNCS__
#define __FUNCS__

// Utility functions
uvec sortNodes (uvec&);
double Correlation (vec&, vec&);
double SignAwareMean (vec&, vec&);

// Network properties
vec WeightedDegree (const mat&, uvec&);
double AverageEdgeWeight (vec&);
vec CorrVector (const mat&, uvec&);
vec SummaryProfile (const mat&, uvec&);
vec NodeContribution (const mat&, uvec&, vec&);
double ModuleCoherence (vec&);
  
#endif // __FUNCS__