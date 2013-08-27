/* This file contains the C++ implementations of the preservation statistics in
 * the corresponding R file, preservation-statistics.R.
 * This file is split into two logical sections: The first part contains 
 * the actual implementations. This is followed by overloaded wrapper functions 
 * whom are needed by Rcpp to generated the actual C++ code (it can't handle the
 * bigmemory dependency, so we hide that with the wrapper functions).
 * See: http://stackoverflow.com/questions/18438291/building-packages-with-rcpp-attributes-not-handled-correctly/
 */


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(bigmemory)]]
//' @useDynLib FastModPres
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

// For debugging
#include <iostream>

// Implementation of meanAdj
NumericVector MeanAdj(XPtr<BigMatrix> pAdjacency, IntegerVector moduleIndices,
                      LogicalVector includeDiagonals) {
  MatrixAccessor<double> adjacency(*pAdjacency);
  NumericVector mean = NumericVector(1);  // Return scalar
  
  // get some useful values
  int moduleSize = moduleIndices.size();
  int totalSize = pAdjacency->nrow()*pAdjacency->ncol();
  
  // Intermediate counters
  int NAcount = 0;
  double total = 0.0; 
  
  // Add to the total sum, handles NAs, and ignores the diagonal if asked.
  for (int i = 0; i < moduleSize; i++) {
    for (int j = 0; j < moduleSize; j++) {
      if ((i != j) || (includeDiagonals[0])) {
        double value = adjacency[moduleIndices[i]-1][moduleIndices[j]-1];
        if (any(is_na(NumericVector(value)))) {
          NAcount += 1;
        } else {
           total += value;
        }
      }
    }
  }
  
  // Divide the total by the number of entries counted.
  if (includeDiagonals[0]) {
    mean = total / (totalSize - NAcount);
  } else {
    mean = total / (totalSize - NAcount - moduleSize);
  }
  
  return mean;
}

// [[Rcpp::export]]
NumericVector MeanAdj(SEXP pAdjacency, IntegerVector moduleIndices, 
                      LogicalVector includeDiagonals) {
  return MeanAdj(XPtr<BigMatrix>(pAdjacency), moduleIndices, includeDiagonals);                  
}
