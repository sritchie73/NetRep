#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
#include <bigmemory/MatrixAccessor.hpp>

/* Implementation of CheckFinite
 *
 * @param xpDat External Pointer for the big.matrix.
 * @param dat MatrixAccessor for the big.matrix. 
 * @return a LogicalVector of length one
 */
template <typename T>
void CheckFinite(XPtr<BigMatrix> xpDat, MatrixAccessor<T> dat) {
  T value;
  int type = xpDat->matrix_type();

  for (unsigned int jj = 0; jj < xpDat->ncol(); jj++) {
    for (unsigned int ii = 0; ii < xpDat->nrow(); ii++) {
      value = dat[jj][ii];
      if (type == 1 && (value == NA_CHAR)) {
        throw Rcpp::exception("'big.matrix' has non-finite values!");
      } else if (type == 2 && (value == NA_SHORT)) {
        throw Rcpp::exception("'big.matrix' has non-finite values!");
      } else if (type == 4 && (value == NA_INTEGER)) {
        throw Rcpp::exception("'big.matrix' has non-finite values!");
      } else if (type == 8 && !R_FINITE(value)) {
        throw Rcpp::exception("'big.matrix' has non-finite values!");
      }
    }
  }
}

//' Check the elements of a `big.matrix`
//' 
//' Are all the values finite? 
//' 
//' @param pDat SEXP container for the pointer to the 
//'   \code{\link[bigmemory]{big.matrix}} to be checked.
//'
//' @rdname chekcFinite-cpp
// [[Rcpp::export]]
void CheckFinite(SEXP pDat) {
  XPtr<BigMatrix> xpDat(pDat);
  
  unsigned short type = xpDat->matrix_type();
  if (type == 1) {
    return CheckFinite(xpDat, MatrixAccessor<char>(*xpDat));
  } else if (type == 2) {
    return CheckFinite(xpDat, MatrixAccessor<short>(*xpDat));
  } else if (type == 4) {
    return CheckFinite(xpDat, MatrixAccessor<int>(*xpDat));
  } else if (type == 8) {
    return CheckFinite(xpDat, MatrixAccessor<double>(*xpDat));
  } else {
    /* We should never get here, unless the underlying implementation of 
    bigmemory changes */
    throw Rcpp::exception("Undefined type for provided big.matrix");
  }
} 
