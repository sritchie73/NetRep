#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
#include <bigmemory/BigMatrix.h>

//' Network subset eigenvector and proportion of variance explained in C++
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
      is_true(any(subsetIndices > xpDat->nrow()))) {
    throw std::out_of_range("Some of requested indices are outside of range!");
  }
  
  // Dispatch function for all types of big.matrix.
  if (xpDat->matrix_type() == 8) {
    // Cast the BigMatrix to an arma::Mat<double>
    mat aDat((double *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false);
    mat U, V;
    vec S;
    uvec subsetRows = as<uvec>(subsetIndices) - 1;
    
    // Get the summary profile for the network subset from the SVD.
    svd_econ(U, S, V, aDat.rows(subsetRows), "right", "dc");
    vec summary(V.col(1));
    
    // Make sure the orientation of the eigenvector matches the orientation of
    // the data
    vec meanExpr(mean(aDat.rows(subsetRows), 1));
    vec sdExpr(stddev(aDat.rows(subsetRows), 0, 1));

    // Scale data in place
    for (unsigned int jj = 0; jj < subsetIndices.size(); jj++) {
      for (unsigned int ii = 0; ii < subsetIndices.size(); ii++) {
        aDat(subsetRows(ii), jj) -= meanExpr(ii);
        aDat(subsetRows(ii), jj) /= sdExpr(ii);
      }   
    }   

    mat ap = cor(mean(aDat.rows(subsetRows)), summary);
    if (ap(0,0) < 0) {                                                                                                                                                                                                                        
      for (unsigned int jj = 0; jj < xpDat->ncol(); jj++) {
        summary(jj) *= -1; 
      }   
    }   

    // Unscale data
    for (unsigned int jj = 0; jj < subsetIndices.size(); jj++) {
      for (unsigned int ii = 0; ii < subsetIndices.size(); ii++) {
        aDat(subsetRows(ii), jj) *= sdExpr(ii);
        aDat(subsetRows(ii), jj) += meanExpr(ii);
      }   
    }   
    
    // The proportion of variance explained is the sum of the squared 
    // correlation between the network subset summary profile, and each of the 
    // variables in the data that correspond to nodes in the network subset.
    mat p = cor(summary, aDat.rows(subsetRows).t());
    for (unsigned int ii = 0; ii < subsetIndices.size(); ii++) {
      p(ii) *= p(ii);
    }
    vec pve(mean(p, 1));
    
    return List::create(
        Named("summaryProfile") = NumericVector(summary.begin(), summary.end())
        Named("propVarExpl") = NumericVector(pve.begin(), pve.end())
      );
  } else {
    throw Rcpp::exception(
      "SVD can only be calculated on a big.matrix whose underlying type is"
      "'double'."
    );
  }
}