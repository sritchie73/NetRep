#define ARMA_USE_LAPACK
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
#include <bigmemory/BigMatrix.h>

/* Implementation of AdjProps
 *
 * @param adj the armadillo compatible adjacency matrix
 * @param subsetIndices indices of the network subset of interest.
 * @return
 *    A List containing:
 *     - The weighted within-subset degree for each node (kIM).
 *     - The mean absolute edge weight of the network subset (meanAdj).
 */
template <typename T>
List AdjProps(const Mat<T>& adj, IntegerVector subsetIndices) {
  // Convert indices to C++ indexing and a class Armadillo can work with.
  // Indices are sorted, because sequential memory access is faster!.
  uvec nodeIdx = sort(as<uvec>(subsetIndices) - 1);


  // We do not want a negative weight to cancel out a positive one, so we take
  // the absolute value.
  Mat<T> dg = diagvec(adj);
  Mat<T> colSums = sum(abs(adj(nodeIdx, nodeIdx))) - abs(dg(nodeIdx)).t();
  
  int n = subsetIndices.size();
  
  // To make sure the resulting KIM vector is in the correct order,
  // order the results to match the original ordering of subsetIndices.
  Function rank("rank"); // Rank only works on R objects like IntegerVector.
  uvec idxRank = as<uvec>(rank(subsetIndices)) - 1;

  Col<T> oKIM = colSums(idxRank);

  return List::create(
    Named("kIM") = oKIM,
    Named("meanAdj") = sum(colSums, 1) / (n*n - n)
  );
}                                                                                                                                                                                                                                          

//' Calculate Mean Adjacency and Intramodular Connectivity
//'
//' @param pAdjacency SEXP container for the pointer to the adjacency matrix.
//' @param subsetIndices indices of the network subset of interest.
//'   
//' @return
//'   A List containing:
//'   \enumerate{
//'     \item{\emph{kIM}:}{The weighted within-subset degree for each node.}
//'     \item{\emph{meanAdj}:}{The mean absolute edge weight of the network subset.}
//'   }
//' @rdname AdjProps-cpp
// [[Rcpp::export]]
List AdjProps(SEXP pAdjacency, IntegerVector subsetIndices) {
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
    return AdjProps(
      arma::Mat<char>((char *)xpAdj->matrix(), xpAdj->nrow(), xpAdj->ncol(), false),
      subsetIndices
    );
  } else if (type == 2) {
    return AdjProps(
      arma::Mat<short>((short *)xpAdj->matrix(), xpAdj->nrow(), xpAdj->ncol(), false),
      subsetIndices
    );
  } else if (type == 4) {
    return AdjProps(
      arma::Mat<int>((int *)xpAdj->matrix(), xpAdj->nrow(), xpAdj->ncol(), false),
      subsetIndices
    );
  } else if (type == 8) {
    return AdjProps(
      arma::Mat<double>((double *)xpAdj->matrix(), xpAdj->nrow(), xpAdj->ncol(), false),
      subsetIndices
    );
  } else {
    /* We should never get here, unless the underlying implementation of
    bigmemory changes */
    throw Rcpp::exception("Undefined type for provided big.matrix");
  }
}
