#ifndef __FUNCS__
#define __FUNCS__

// [[Rcpp::depends(BH, RcppArmadillo)]]
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