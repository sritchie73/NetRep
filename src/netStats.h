#ifndef __FUNCS__
#define __FUNCS__

#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#define ARMA_NO_DEBUG
#define ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_USE_CXX11

#include <RcppArmadillo.h>

// Utility functions
arma::uvec SortNodes (unsigned int *, unsigned int);
double Correlation (double *, double *, unsigned int);
double SignAwareMean (double *, double *, unsigned int);

// Network properties
arma::vec WeightedDegree (double *, unsigned int, unsigned int *, unsigned int);
double AverageEdgeWeight (double *, unsigned int);
arma::vec CorrVector (double *, unsigned int, unsigned int *, unsigned int);
arma::vec SummaryProfile (double *, unsigned int, unsigned int, unsigned int *, unsigned int);
arma::vec NodeContribution (double *, unsigned int, unsigned int, unsigned int *, unsigned int, double *);
double ModuleCoherence (double *, unsigned int);
  
#endif // __FUNCS__