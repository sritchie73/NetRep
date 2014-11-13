#define ARMA_USE_LAPACK
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
#include <bigmemory/BigMatrix.h>

template <typename T>
List RangeSubset(const Mat<T>& dat, IntegerVector subsetIndices) {
  // Convert indices to C++ indexing and a class Armadillo can work with.
  // Indices are sorted, because sequential memory access is faster!.
  uvec nodeIdx = sort(as<uvec>(subsetIndices) - 1);
  
  Mat<T> datMax = max(max(dat.cols(nodeIdx), 0), 1);
  Mat<T> datMin = min(min(dat.cols(nodeIdx), 0), 1);
  
  return List::create(
    Named("min") = datMin,
    Named("max") = datMax
  );
}

template <typename T>
List BigRange(const Mat<T>& dat) {
  Mat<T> datMax = max(max(dat, 0), 1);
  Mat<T> datMin = min(min(dat, 0), 1);
  
  return List::create(
    Named("min") = datMin,
    Named("max") = datMax
  );
}

//' Get the range of a big.matrix
//' 
//' @description
//'  \code{RangeSubset}: get the range of values in the column-subset of a 
//'  big.matrix.
//' 
//' @param pDat SEXP container for the pointer to the data matrix to be scaled.
//' @param subsetIndices indices of the network subset of interest in 
//'   \code{pDat}.
//'   
//' @rdname range-cpp
// [[Rcpp::export]]
List RangeSubset(SEXP pDat, IntegerVector subsetIndices) {
  XPtr<BigMatrix> xpDat(pDat);
  
  // Make sure we're not indexing out of range.
  if (is_true(any(subsetIndices <= 0)) ||
      is_true(any(subsetIndices > xpDat->ncol()))) {
    throw std::out_of_range("requested index outside of range!");
  }

  //  Dispatch function for all types of big.matrix.
  unsigned short type = xpDat->matrix_type();
  if (type == 1) {
    return RangeSubset(
      arma::Mat<char>((char *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false),
      subsetIndices
    );
  } else if (type == 2) {
    return RangeSubset(
      arma::Mat<short>((short *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false),
      subsetIndices
    );
  } else if (type == 4) {
    return RangeSubset(
      arma::Mat<int>((int *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false),
      subsetIndices
    );
  } else if (type == 8) {
    return RangeSubset(
      arma::Mat<double>((double *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false),
      subsetIndices
    );
  } else {
    /* We should never get here, unless the underlying implementation of
    bigmemory changes */
    throw Rcpp::exception("Undefined type for provided big.matrix");
  }
} 

//' @name range-cpp
//' @description
//'   \code{BigRange}: get the range of values in a big.matrix
//' 
// [[Rcpp::export]]
List BigRange(SEXP pDat) {
  XPtr<BigMatrix> xpDat(pDat);

  //  Dispatch function for all types of big.matrix.
  unsigned short type = xpDat->matrix_type();
  if (type == 1) {
    return BigRange(
      arma::Mat<char>((char *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false)
    );
  } else if (type == 2) {
    return BigRange(
      arma::Mat<short>((short *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false)
    );
  } else if (type == 4) {
    return BigRange(
      arma::Mat<int>((int *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false)
    );
  } else if (type == 8) {
    return BigRange(
      arma::Mat<double>((double *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false)
    );
  } else {
    /* We should never get here, unless the underlying implementation of
    bigmemory changes */
    throw Rcpp::exception("Undefined type for provided big.matrix");
  }
} 

