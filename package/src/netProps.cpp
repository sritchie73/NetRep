#define ARMA_USE_LAPACK
#define ARMA_NO_DEBUG
#define ARMA_DONT_USE_CXX11

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
#include <bigmemory/BigMatrix.h>

/* Implementation of NetProps
 *
 * @param adj the armadillo compatible adjacency matrix
 * @param subsetIndices indices of the network subset to compute the mean
 *   adjacency on
 * @return
 *    A List containing:
 *     - The weighted within-subset degree for each node (kIM).
 *     - The maximum adjacency ratio for each node (MAR).
 *     - The mean absolute edge weight of the network subset (meanAdj).
 *     - The mean within-subset degree (meanKIM).
 *     - The mean maximum adjacency ratio (meanMAR).
 */
template <typename T>
List NetProps(const Mat<T>& adj, IntegerVector subsetIndices) {
  // Convert indices to C++ indexing and a class Armadillo can work with
  uvec nodeIdx = as<uvec>(subsetIndices) - 1;

  // We do not want a negative weight to cancel out a positive one, so we take
  // the absolute value.
  Mat<T> dg = diagvec(adj);
  Mat<T> colSums = sum(abs(adj(nodeIdx, nodeIdx))) - abs(dg(nodeIdx)).t();
  Mat<T> meanKIM = mean(colSums, 1); // This will be length 1

  Mat<T> sqSums = sum(square(adj(nodeIdx, nodeIdx))) - square(dg(nodeIdx)).t();
  Mat<T> MAR = sqSums / colSums;
  Mat<T> meanMAR = mean(MAR, 1); // This will be length 1
  
  int n = subsetIndices.size();

  return List::create(
    Named("kIM") = colSums,
    Named("MAR") = MAR,
    Named("meanAdj") = sum(colSums, 1) / (n*n - n), 
    Named("meanKIM") = meanKIM,                                                                                                                                                                                                            
    Named("meanMAR") = meanMAR
  );
}                                                                                                                                                                                                                                          

//' Calculate Network Properties
//'
//' @param pAdjacency SEXP container for the pointer to the adjacency matrix
//' @param subsetIndices indices of the subset of the network to calculate
//'   the mean adjacency for.
//'   
//' @return
//'   A List containing:
//'   \enumerate{
//'     \item{\emph{kIM}:}{The weighted within-subset degree for each node.}
//'     \item{\emph{MAR}:}{The maximum adjacency ratio for each node.}
//'     \item{\emph{meanAdj}:}{The mean absolute edge weight of the network subset.}
//'     \item{\emph{meanKIM}:}{The mean within-subset degree.}
//'     \item{\emph{meanMAR}:}{The mean maximum adjacency ratio.}
//'   }
//' @rdname netProps-cpp
// [[Rcpp::export]]
List NetProps(SEXP pAdjacency, IntegerVector subsetIndices) {
  XPtr<BigMatrix> xpAdj(pAdjacency);

  // Make sure we're not indexing out of range.
  if (is_true(any(subsetIndices <= 0)) ||
      is_true(any(subsetIndices > xpAdj->ncol())) ||
      is_true(any(subsetIndices > xpAdj->nrow()))) {
    throw std::out_of_range("Requested index outside of range!");
  }

  //  Dispatch function for all types of big.matrix.
  unsigned short type = xpAdj->matrix_type();
  if (type == 1) {
    return NetProps(
      arma::Mat<char>((char *)xpAdj->matrix(), xpAdj->nrow(), xpAdj->ncol(), false),
      subsetIndices
    );
  } else if (type == 2) {
    return NetProps(
      arma::Mat<short>((short *)xpAdj->matrix(), xpAdj->nrow(), xpAdj->ncol(), false),
      subsetIndices
    );
  } else if (type == 4) {
    return NetProps(
      arma::Mat<int>((int *)xpAdj->matrix(), xpAdj->nrow(), xpAdj->ncol(), false),
      subsetIndices
    );
  } else if (type == 8) {
    return NetProps(
      arma::Mat<double>((double *)xpAdj->matrix(), xpAdj->nrow(), xpAdj->ncol(), false),
      subsetIndices
    );
  } else {
    /* We should never get here, unless the underlying implementation of
    bigmemory changes */
    throw Rcpp::exception("Undefined type for provided big.matrix");
  }
}
