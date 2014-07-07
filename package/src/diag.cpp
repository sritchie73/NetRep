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

//' big.matrix diagonals
//' 
//' @param pBigMat a \code{numeric} \code{\link[bigmemory]{big.matrix}}
//' @return a vector containing the diagonal of \code{pBigMat}
// [[Rcpp::export]]
NumericVector GetDiag(SEXP pBigMat) {
  return GetDiag(XPtr<BigMatrix>(pBigMat));                  
}

//' big.matrix diagonals
//' 
//' @param pBigMat a \code{numeric} \code{\link[bigmemory]{big.matrix}}
//' @param value either a single \code{numeric} value or \code{numeric} vector 
//'   of length equal to the current diagonal.
// [[Rcpp::export]]
void SetDiag(SEXP pBigMat, NumericVector value) {
  return SetDiag(XPtr<BigMatrix>(pBigMat), value);                  
}
