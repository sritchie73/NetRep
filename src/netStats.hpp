#ifndef __FUNCS__
#define __FUNCS__

#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#define ARMA_NO_DEBUG
#define ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_USE_CXX11

#include <RcppArmadillo.h>

// Utility functions
arma::mat Scale (const arma::mat&);
arma::uvec sortNodes (arma::uvec&);
double Correlation (arma::vec&, arma::vec&);
double SignAwareMean (arma::vec&, arma::vec&);

// Network properties
arma::vec WeightedDegree (const arma::mat&, arma::uvec&);
double AverageEdgeWeight (arma::vec&);
arma::vec CorrVector (const arma::mat&, arma::uvec&);
arma::vec SummaryProfile (const arma::mat&, arma::uvec&);
arma::vec NodeContribution (const arma::mat&, arma::uvec&, arma::vec&);
double ModuleCoherence (arma::vec&);
  
#endif // __FUNCS__