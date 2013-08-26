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
// TODO: Handle NA's in entries.
NumericVector MeanAdj(XPtr<BigMatrix> pAdjacency, IntegerVector moduleIndices,
                      LogicalVector includeDiagonals) {
  MatrixAccessor<double> adjacency(*pAdjacency);
  int moduleSize = moduleIndices.size();
  int totalSize = pAdjacency->nrow()*pAdjacency->ncol();
  NumericVector mean = NumericVector(1);
  
  double total = 0.0; 
  
  for (int i = 0; i < moduleSize; i++) {
    for (int j = 0; j < moduleSize; j++) {
      if ((i != j) || (includeDiagonals[0])) {
        total += adjacency[moduleIndices[i]-1][moduleIndices[j]-1];
      }
    }
  }
  
  if (includeDiagonals[0]) {
    mean = total / totalSize;
  } else {
    mean = total / (totalSize - moduleSize);
  }
  
  return mean;
}

// [[Rcpp::export]]
NumericVector MeanAdj(SEXP pAdjacency, IntegerVector moduleIndices, 
                      LogicalVector includeDiagonals) {
  return MeanAdj(XPtr<BigMatrix>(pAdjacency), moduleIndices, includeDiagonals);                  
}
