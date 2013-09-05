/** 
 * This file contains the C++ implementation of the kIM calculations on a 
 * big.matrix pointer.
 */

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

// This is colSums over a subset of a matrix.
NumericVector KIM(XPtr<BigMatrix> pAdjacency, IntegerVector moduleIndices) {
  
  MatrixAccessor<double> adjacency(*pAdjacency);
  NumericVector kIM = NumericVector(moduleIndices.size());  // Return scalar
  
  // temporary value holder
  NumericVector value = NumericVector(1); 
  
  // Get some useful values
  int ncol = pAdjacency->ncol();
  int nrow = pAdjacency->nrow();
  
  // Make sure we're not indexing out of range.
  if (is_true(any(moduleIndices <= 0)) || is_true(any(moduleIndices > ncol)) ||
      is_true(any(moduleIndices > nrow))) {
    throw std::out_of_range("Requested index outside of range!");
  }

  for (int i = 0; i < moduleIndices.size(); i++) {
    kIM[i] = 0;
    for (int j = 0; j < moduleIndices.size(); j++) {
      value = adjacency[moduleIndices[i]-1][moduleIndices[j]-1];
      if (all(!is_na(value))) {
        kIM[i] += value[0];
      }
    }
  }
  
  return kIM;
}

// [[Rcpp::export]]
NumericVector KIM(SEXP pAdjacency, IntegerVector moduleIndices) {
  return KIM(XPtr<BigMatrix>(pAdjacency), moduleIndices);                  
}
