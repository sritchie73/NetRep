/** 
 * This file contains the C++ implementation of the kIM calculations on a 
 * big.matrix pointer.
 */

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

#include "utilities.h"

// This is colSums over a subset of a matrix.
NumericVector KIM(XPtr<BigMatrix> pAdjacency, IntegerVector moduleIndices,
                  LogicalVector allNodes) {
  
  MatrixAccessor<double> adjacency(*pAdjacency);
  NumericVector kIM = NumericVector(moduleIndices.size());  // Return scalar
  
  // temporary value holder
  NumericVector value = NumericVector(1); 

  for (int i = 0; i < moduleIndices.size(); i++) {
    kIM[i] = 0;
    // Get connectivity to all nodes? or just module nodes?
    if(allNodes[0]) {
      for (int j = 0; j < pAdjacency->nrow(); j++) {
        value = safeAccessor(pAdjacency, j, moduleIndices[i]-1);
        if (all(!is_na(value))) {
          kIM[i] += value[0];
        }
      }
    } else {
      for (int j = 0; j < moduleIndices.size(); j++) {
        value = safeAccessor(pAdjacency, moduleIndices[j]-1, moduleIndices[i]-1);
        if (all(!is_na(value))) {
          kIM[i] += value[0];
        }
      }
    }
  }
  
  return kIM;
}

// [[Rcpp::export]]
NumericVector KIM(SEXP pAdjacency, IntegerVector moduleIndices, 
                  LogicalVector allNodes) {
  return KIM(XPtr<BigMatrix>(pAdjacency), moduleIndices, allNodes);                  
}
