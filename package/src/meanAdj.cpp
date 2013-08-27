/** 
 * This file contains the C++ implementation of the meanAdj statistic on a 
 * big.matrix pointer.
 */

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(bigmemory)]]
//' @useDynLib FastModPres
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

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
  NumericVector value = NumericVector(1);
  
  // Add to the total sum, handles NAs, and ignores the diagonal if asked.
  for (int i = 0; i < moduleSize; i++) {
    for (int j = 0; j < moduleSize; j++) {
      if ((i != j) || (includeDiagonals[0])) {
        value[0] = adjacency[moduleIndices[i]-1][moduleIndices[j]-1];
        if (any(is_na(value))) {
          NAcount += 1;
        } else { 
          total += value[0];
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
