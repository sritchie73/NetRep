#define ARMA_USE_LAPACK
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
#include <bigmemory/BigMatrix.h>

//' Scale a matrix by its rows
//' 
//' @param pDat SEXP container for the pointer to the data matrix to be scaled.
//' @param spDat SEXP container for the pointer to the pre-initialised
//'   \code{\link[bigmemory]{big.matrix}} that the scaled version of \code{pDat}
//'   will be stored in.
//' @rdname scale-cpp
// [[Rcpp::export]]
void Scale(SEXP pDat, SEXP spDat) {
  XPtr<BigMatrix> xpDat(pDat);
  XPtr<BigMatrix> xspDat(spDat);
  
  // Make sure both matrices are concordant.
  if (xpDat->ncol() != xspDat->ncol() ||
      xpDat->nrow() != xspDat->nrow()) {
    throw Rcpp::exception(
        "The results matrix must have the same dimensions as the data matrix!"
      );
  }
  if (xpDat->matrix_type() != xspDat->matrix_type()) {
    throw Rcpp::exception(
        "The results matrix must have the same 'type' as the data matrix."
      );
  }
  
  if (xpDat->matrix_type() == 8) {
    // Cast both matrices to Armadillo types
    mat aDat((double *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false);
    mat asDat((double *)xspDat->matrix(), xspDat->nrow(), xspDat->ncol(), false);
    
    // Get the mean and std for each column (network node)
    rowvec meanExpr(mean(aDat));
    rowvec sdExpr(stddev(aDat));
  
    // Store scaled data
    for (unsigned int jj = 0; jj < xpDat->ncol(); jj++) {
      asDat.col(jj) = (aDat.col(jj) - meanExpr(jj))/sdExpr(jj);
    }
  } else {
    throw Rcpp::exception(
      "Network statistics calculated from the underlying data only work on"
      " big.matrix objects of type 'double'."
    );
  }
} 
