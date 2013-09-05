/** 
 * This file contains the C++ implementation of the meanAdj statistic on a 
 * big.matrix pointer.
 */
 
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

//' @useDynLib FastModPres


NumericVector MeanAdj(XPtr<BigMatrix> pAdjacency, IntegerVector moduleIndices) {
  NumericVector mean = NumericVector(1);  // Return scalar
  MatrixAccessor<double> adjacency(*pAdjacency);
  
  // get some useful values
  int moduleSize = moduleIndices.size();
  int ncol = pAdjacency->ncol();
  int nrow = pAdjacency->nrow();
  
  // Intermediate counters
  int NAcount = 0;
  double total = 0.0;
  double value;
  
  // Make sure we're not indexing out of range.
  if (is_true(any(moduleIndices <= 0)) || is_true(any(moduleIndices > nrow)) ||
      is_true(any(moduleIndices > ncol))) {
    throw std::out_of_range("Requested index outside of range!");
  }
  
  // Add to the total sum while handling NAs
  for (int i = 0; i < moduleSize; i++) {
    for (int j = 0; j < moduleSize; j++) {
     value = adjacency[moduleIndices[i]-1][moduleIndices[j]-1];
      if (R_IsNA(value)) {
        NAcount += 1;
      } else { 
        total += value;
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
