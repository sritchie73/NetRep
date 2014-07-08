// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(BH, bigmemory)]]
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

/* Generic implementation of GetDiag
 *
 * @param pMat The pointer to the big.matrix object
 * @param mat The MatrixAccessor object, allowing us to access the values
 *  within the big.matrix object.
 *
 * @return The vector of values stored in the diagonal of \code{matVals}
 */
template <typename T>
NumericVector GetDiag(XPtr<BigMatrix> pMat, MatrixAccessor<T> mat) {
  int diagLength = GetDiagLength(pMat);
  NumericVector resultsVector(diagLength);
  for (int i = 0; i < diagLength; i++) {
      resultsVector[i] = mat[i][i];
  }
  return resultsVector;
}

//' big.matrix diagonals
//' 
//' @param pBigMat a \code{numeric} \code{\link[bigmemory]{big.matrix}}
//' @return a vector containing the diagonal of \code{pBigMat}
// [[Rcpp::export]]
NumericVector GetDiag(SEXP pBigMat) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
    case 1:
      return GetDiag(xpMat, MatrixAccessor<char>(*xpMat));
    case 2:
      return GetDiag(xpMat, MatrixAccessor<short>(*xpMat));
    case 4:
      return GetDiag(xpMat, MatrixAccessor<int>(*xpMat));
    case 8:
      return GetDiag(xpMat, MatrixAccessor<double>(*xpMat));
    default:
      /* We should never get here, unless the underlying implementation of 
         bigmemory changes */
      throw Rcpp::exception("Undefined type for provided big.matrix");
  }
}

/* Generic implementation of SetDiag
 *
 * @param pMat The pointer to the big.matrix object
 * @param mat The MatrixAccessor object, allowing us to access the values
 *  within the big.matrix object.
 * @param value The NumericVector
 *
 * @return The vector of values stored in the diagonal of \code{matVals}
 */
template <typename T>
void SetDiag(XPtr<BigMatrix> pMat, MatrixAccessor<T> mat, NumericVector value) {
  int diagLength = GetDiagLength(pMat);
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
//' @param value either a single \code{numeric} value or \code{numeric} vector 
//'   of length equal to the current diagonal.
// [[Rcpp::export]]
void SetDiag(SEXP pBigMat, NumericVector value) {
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
    case 1:
      return SetDiag(xpMat, MatrixAccessor<char>(*xpMat), value);
    case 2:
      return SetDiag(xpMat, MatrixAccessor<short>(*xpMat), value);
    case 4:
      return SetDiag(xpMat, MatrixAccessor<int>(*xpMat), value);
    case 8:
      return SetDiag(xpMat, MatrixAccessor<double>(*xpMat), value);
    default:
      /* We should never get here, unless the underlying implementation of 
         bigmemory changes */
      throw Rcpp::exception("Undefined type for provided big.matrix");
  }              
}
