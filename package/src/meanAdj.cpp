/** 
 * This file contains the C++ implementation of the meanAdj statistic on a 
 * big.matrix pointer.
 */
 
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

#include "utilities.h"

NumericVector MeanAdj(XPtr<BigMatrix> pAdjacency, IntegerVector moduleIndices) {
  NumericVector mean = NumericVector(1);  // Return scalar
  
  // get some useful values
  int moduleSize = moduleIndices.size();
  
  // Intermediate counters
  int NAcount = 0;
  double total = 0.0;
  NumericVector value = NumericVector(1);
  
  // Add to the total sum, handles NAs, and ignores the diagonal if asked.
  for (int i = 0; i < moduleSize; i++) {
    for (int j = 0; j < moduleSize; j++) {
      value = safeAccessor(pAdjacency, moduleIndices[i]-1, moduleIndices[j]-1);
      if (any(is_na(value))) {
        NAcount += 1;
      } else { 
        total += value[0];
      }
    }
  }

  mean = total / (moduleSize * moduleSize - NAcount);
  
  return mean;
}

// [[Rcpp::export]]
NumericVector MeanAdj(SEXP pAdjacency, IntegerVector moduleIndices) {
  return MeanAdj(XPtr<BigMatrix>(pAdjacency), moduleIndices);                  
}
