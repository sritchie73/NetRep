#include <Rcpp.h>
using namespace Rcpp;

//' Network subset eigenvector and proportion of variance explained in C++
//' 
//' @param pAdjacency SEXP container for the pointer to the adjacency matrix
//' @param pDat SEXP container for the pointer to the data matrix that 
//'   \code{pAdjacency} was constructed from.
//' @param adjIndices indices of the network subset of interest in
//'   \code{pAdjacency}.
//' @param datIndices indices of the network subset of interest in \code{pDat}.
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
// [[Rcpp::export]]
NumericVector SvdProps(
  SEXP pAdjacency, SEXP pDat, IntegerVector adjIndices, IntegerVector datIndices
) {
  XPtr<BigMatrix> xpAdj(pAdjacency);
  XPtr<BigMatrix> xpDat(pDat);
  
  // Make sure we're not indexing out of range.
  if (is_true(any(adjIndices <= 0)) || 
      is_true(any(adjIndices > xpAdj->ncol())) ||
      is_true(any(adjIndices > xpAdj->nrow()))) {
    throw std::out_of_range(
      "Some of the requested indices for network subset are outside of the "
      "given adjacency matrix."
    );
  }
  if (is_true(any(datIndices <= 0)) || 
      is_true(any(datIndices > xpDat->ncol()))) {
    throw std::out_of_range(
      "Some of the requested indices for network subset are outside of the "
      "given data matrix."
    );
  }
  
  
}