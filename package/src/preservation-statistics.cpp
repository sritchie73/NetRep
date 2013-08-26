
#include <Rcpp.h>

// [[Rcpp::depends(bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

Rcpp::NumericVector BigColSums(Rcpp::XPtr<BigMatrix> pBigMat) {

    // Create the matrix accessor so we can get at the elements of the matrix.
    MatrixAccessor<double> ma(*pBigMat);
  
    // Create the vector we'll store the column sums in.
    Rcpp::NumericVector colSums(pBigMat->ncol());
    for (size_t i=0; i < pBigMat->ncol(); ++i)
        colSums[i] = std::accumulate(ma[i], ma[i]+pBigMat->nrow(), 0.0);
    return colSums;
}

// Wrapper to be exported
// [[Rcpp::export]]
Rcpp::NumericVector BigColSums( SEXP pBigMat ) {
  return BigColSums( Rcpp::XPtr<BigMatrix>(pBigMat) );
}