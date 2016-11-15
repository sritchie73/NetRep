#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#define ARMA_NO_DEBUG
#define ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_USE_CXX11

#include <RcppArmadillo.h>

///' Check whether there are any non-finite values in a matrix
///'
///' The C++ functions will not work with NA values, and the calculation of the
///' summary profile will take a long time to run before crashing.
///'
///' @param matPtr matrix to check.
///' 
///' @return
///'  Throws an error if any \code{NA}, \code{NaN}, \code{Inf}, or \code{-Inf}
///'  values are found, otherwise returns silently.
///' 
///' @keywords internal
// [[Rcpp::export]]
void CheckFinite(Rcpp::NumericMatrix matPtr) {
  arma::mat mat = arma::mat(matPtr.begin(), matPtr.nrow(), matPtr.ncol(), false, true);
  arma::uvec nonFiniteIdx = arma::find_nonfinite(mat);
  if (nonFiniteIdx.n_elem > 0) {
    throw Rcpp::exception("matrices cannot have non-finite or missing values");
  } 
}
