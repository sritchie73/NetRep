#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
#include <bigmemory/BigMatrix.h>


/* Implementation for SvdProps
 *
 * @param dat Matrix Accessor for the data matrix.
 * @param subsetIndices indices of the network subset in 'dat'.
 * @return
 *  A List containing two 'NumericVector's (the first eigenvector of the network 
 *  subset in 'dat', and the proportion of variance in 'dat' it explains).
 */
template <typename T>
List SvdProps(const arma::Mat<T>& aDat, IntegerVector subsetIndices) {
  return List::create(CharacterVector("Success"));
} 

//' Network subset eigenvector and proportion of variance explained in C++
//' 
//' C++ Dispatch Function
//' 
//' @param pDat SEXP container for the pointer to the data matrix used in 
//'   network construction.
//' @param subsetIndices indices of the network subset of interest in 
//'   \code{pDat}.
//' 
//' @return
//'  A list whose first element is the first eigenvector of the singular value
//'  decomposition for the network subset, and whose second element is the 
//'  proportion of the variance in the corresponding subset of \code{pDat} the 
//'  eigenvector explains.
//'  
//' @details
//'  The sign of the returned eigenvector is modified to match the average of
//'  \code{pDat}. This is to match the behaviour of 
//'  \code{\link[WGCNA]{moduleEigengenes}} in the \code{\link{WGCNA}} package.
//'  
//'  The two returned properties are bundled together into one function because
//'  the calculation of the proportion of variance requires much of the same
//'  underlying intermediate calculations that obtaining the first eigenvector
//'  requires.
//' 
//' @import RcppArmadillo
//'  
// [[Rcpp::export]]
List SvdProps(
  SEXP pDat, IntegerVector subsetIndices
) {
  XPtr<BigMatrix> xpDat(pDat);
  
  // Make sure we're not indexing out of range.
  if (is_true(any(subsetIndices <= 0)) || 
      is_true(any(subsetIndices > xpDat->ncol()))) {
    throw std::out_of_range(
      "Some of the requested indices for network subset are outside of the "
      "given data matrix."
    );
  }
  
  // Dispatch function for all types of big.matrix.
  unsigned short datType = xpDat->matrix_type();
  if (datType == 1) {
    return SvdProps(
      Mat<char>((char *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false),
      subsetIndices
    );
  } else if (datType == 2) {
    return SvdProps(
      Mat<short>((short *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false),
      subsetIndices
    );  
  } else if (datType == 4) {
    return SvdProps(
      Mat<int>((int *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false),
      subsetIndices
    );
  } else if (datType == 8) {
    return SvdProps(
      Mat<double>((double *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false),
      subsetIndices
    );
  } else {
    /* We should never get here, unless the underlying implementation of 
    bigmemory changes */
    throw Rcpp::exception("Undefined type for provided data big.matrix");
  }
}