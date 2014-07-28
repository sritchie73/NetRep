#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
#include <bigmemory/MatrixAccessor.hpp>

/* Implementation of AllFinite
 *
 * @param xpDat External Pointer for the big.matrix.
 * @param dat MatrixAccessor for the big.matrix. 
 * @return a LogicalVector of length one
 */
template <typename T>
LogicalVector AllFinite(XPtr<BigMatrix> xpDat, MatrixAccessor<T> dat) {
  T value;
  int type = xpDat->matrix_type();

  for (unsigned int jj = 0; jj < xpDat->ncol(); jj++) {
    for (unsigned int ii = 0; ii < xpDat->nrow(); ii++) {
      value = dat[jj][ii];
      if (type == 1 && !(value == NA_CHAR)) {
        return wrap(false);
      } else if (type == 2 && !(value == NA_SHORT)) {
        return wrap(false);
      } else if (type == 4 && !(value == NA_INTEGER)) {
        return wrap(false);
      } else if (type == 8 && !R_FINITE(value)) {
        return wrap(false);
      }
    }
  }
  return wrap(true);
}

//' Check the elements of a `big.matrix`
//' 
//' Are all the values finite? 
//' 
//' @param pDat SEXP container for the pointer to the 
//'   \code{\link[bigmemory]{big.matrix}} to be checked.
//' @return
//'   \code{TRUE} if all values are finite, \code{FALSE} otherwise
//' @rdname allFinite-cpp
// [[Rcpp::export]]
LogicalVector AllFinite(SEXP pDat) {
  XPtr<BigMatrix> xpDat(pDat);
  
  unsigned short type = xpDat->matrix_type();
  if (type == 1) {
    return AllFinite(xpDat, MatrixAccessor<char>(*xpDat));
  } else if (type == 2) {
    return AllFinite(xpDat, MatrixAccessor<short>(*xpDat));
  } else if (type == 4) {
    return AllFinite(xpDat, MatrixAccessor<int>(*xpDat));
  } else if (type == 8) {
    return AllFinite(xpDat, MatrixAccessor<double>(*xpDat));
  } else {
    /* We should never get here, unless the underlying implementation of 
    bigmemory changes */
    throw Rcpp::exception("Undefined type for provided big.matrix");
  }
} 
