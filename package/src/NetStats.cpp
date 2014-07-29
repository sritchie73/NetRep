#define ARMA_USE_LAPACK
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
#include <bigmemory/BigMatrix.h>

//' Calculate Network Statistics
//'
//' For some statistics it does not make sense to calculate the necessary
//' components in advance due to large memory overhead, or logic that doesn't
//' separate nicely. This function deals with those statistics.
//'
//' @param pAdjD,pAdjT SEXP containers for the pointers to the
//'  \code{\link[bigmemory]{big.matrix}} objects holding the \emph{discovery}
//'  and \emph{test} networks respectively.
//' @param discIndices,testIndices indices of the network subset of interest in
//'   the \emph{discovery} and \emph{test} networks respectively.
//'
//' @return
//'   A vector containing:
//'   \enumerate{
//'     \item{\emph{cor.adj}:}{
//'       The correlation between the edge weights for both networks.
//'     }
//'   }
//' @rdname netStats-cpp
// [[Rcpp::export]]
List NetStats(
  SEXP pAdjD, IntegerVector discIndices,
  SEXP pAdjT, IntegerVector testIndices
) {
  XPtr<BigMatrix> xpAdjD(pAdjD);
  XPtr<BigMatrix> xpAdjT(pAdjT);

  // Make sure we're not indexing out of range.
  if (is_true(any(discIndices <= 0)) ||
      is_true(any(discIndices > xpAdjD->ncol())) ||
      is_true(any(discIndices > xpAdjD->nrow()))) {
    throw std::out_of_range(
      "Some of the requested indices for the discovery network are out of range!"
    );
  }
  if (is_true(any(testIndices <= 0)) ||
     is_true(any(testIndices > xpAdjT->ncol())) ||
     is_true(any(testIndices > xpAdjT->nrow()))) {
    throw std::out_of_range(
      "Some of the requested indices for the test network are out of range!"
    );
  }

  //  Dispatch function for all types of big.matrix.
  unsigned short typeD = xpAdjD->matrix_type();
  unsigned short typeT = xpAdjT->matrix_type();
  if (typeD != typeT) {
    throw Rcpp::exception(
      "Both big.matrix objects must have the same underlying type!"
    );
  }
  if (typeD == 8) {
    mat discAdj((double *)xpAdjD->matrix(), xpAdjD->nrow(), xpAdjD->ncol(), false);
    mat testAdj((double *)xpAdjT->matrix(), xpAdjT->nrow(), xpAdjT->ncol(), false);
    uvec discIdx = as<uvec>(discIndices) - 1;
    uvec testIdx = as<uvec>(testIndices) - 1;
    mat p = cor(
      vectorise(discAdj(discIdx, discIdx)),
      vectorise(testAdj(testIdx, testIdx))
    );
    return List::create(
      Named("cor.adj") = p
    );
  } else {
    Function warning("warning");
    warning("cor.adj can only be calculated for big.matrix objects of type"
        " \"double\"");
    return List::create(
      Named("cor.adj") = NA_REAL
    );
  }
}
