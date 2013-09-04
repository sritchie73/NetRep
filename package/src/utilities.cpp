/**
 * This file contains utility functions for dealing with big.matrix objects
 */
 
#include <Rcpp.h>
using namespace Rcpp;

//' @useDynLib FastModPres

// [[Rcpp::depends(bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

#include "utilities.h"

NumericVector safeAccessor(XPtr<BigMatrix> pBigMat, int row, int column) {
  MatrixAccessor<double> mat(*pBigMat);
  NumericVector value(1);
  
  if ((row >= 0) && (row < pBigMat->nrow()) && 
      (column >= 0) && (column < pBigMat->ncol())) {
    value = mat[column][row];
    return value;
  } else {
    throw std::out_of_range("Requested index outside of range!");
  }
}
