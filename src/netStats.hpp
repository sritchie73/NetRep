#ifndef __FUNCS__
#define __FUNCS__

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

// Utility functions
mat Scale (const mat&);
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