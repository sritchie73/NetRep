/** 
 * Get and Set the diagonal of a big.matrix.
 */

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

// Get the length of a diagonal of a big.matrix.
int GetDiagLength(XPtr<BigMatrix> pBigMat) {
  if(pBigMat->ncol() < pBigMat->nrow()) {
    return pBigMat->ncol();
  } else {
    return pBigMat->nrow();
  }
}

NumericVector GetDiag(XPtr<BigMatrix> pBigMat) {
  MatrixAccessor<double> mat(*pBigMat);
  int diagLength = GetDiagLength(pBigMat);
  NumericVector diag(diagLength);
  
  for (int i = 0; i < diagLength; i++) {
      diag[i] = mat[i][i];
  }
  
  return diag;
}

void SetDiag(XPtr<BigMatrix> pBigMat, NumericVector value) {
  MatrixAccessor<double> mat(*pBigMat);
  int diagLength = GetDiagLength(pBigMat);
  int valIndex = 0;
  
  if ((value.size() != diagLength) && (value.size() != 1)) {
    throw std::range_error("replacement diagonal has wrong length");
  }
  
  for (int i = 0; i < diagLength; i++) {
    if (value.size() != 1) valIndex = i;
    mat[i][i] = value[valIndex];
  }
}

// [[Rcpp::export]]
NumericVector GetDiag(SEXP pBigMat) {
  return GetDiag(XPtr<BigMatrix>(pBigMat));                  
}

// [[Rcpp::export]]
void SetDiag(SEXP pBigMat, NumericVector value) {
  return SetDiag(XPtr<BigMatrix>(pBigMat), value);                  
}
