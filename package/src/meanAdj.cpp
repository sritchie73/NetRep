#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

#include <cmath>

/* Implementation of Mean Adjacency
 * 
 * @param xpAdj External Pointer for the adjacency matrix.
 * @param adj Matrix Accessor for the adjacency matrix.
 * @param subsetIndices indices of the network subset to compute the mean 
 *   adjacency on
 * @param undirected logical; if TRUE only the lower half of 'adj' is used in
 *   the calculations
 * @return a length one NumericVector
 */
template <typename T>
NumericVector MeanAdj(
  XPtr<BigMatrix> xpAdj, MatrixAccessor<T> adj, IntegerVector subsetIndices,
  LogicalVector undirected
) {
  NumericVector mean = NumericVector(1);  // Return scalar
  
  // get some useful values
  int subsetSize = subsetIndices.size();
  int type = xpAdj->matrix_type();
  
  // Intermediate counters
  int NAcount = 0;
  double total = 0.0;
  double value;
  
  // Add to the total sum while handling NAs
  for (int jj = 0; jj < subsetSize; jj++) {
    unsigned int ii = 0;
    if (undirected[0]) {
      ii = jj+1;
    }
    for (int ii = 0; ii < subsetSize; ii++) {
      value = adj[subsetIndices[jj]-1][subsetIndices[ii]-1];
      if (type == 1 && !(value == NA_CHAR)) {
        total += abs(value);
      } else if (type == 2 && !(value == NA_SHORT)) {
        total += abs(value);
      } else if (type == 4 && !(value == NA_INTEGER)) {
        total += abs(value);
      } else if (type == 8 && !ISNA(value) && !ISNAN(value)) {
        total += abs(value);
      } else {
        NAcount += 1;
      }
    }
  }

  mean = total / (subsetSize * subsetSize - NAcount);
  
  return mean;
}

//' C++ implementation of Mean Adjacency
//' 
//' See the \link[=meanAdj]{wrapper function} for documentation.
//' 
//' @param pAdjacency SEXP container for the pointer to the adjacency matrix
//' @param subsetIndices indices of the subset of the network to calculate
//'   the mean adjacency for.
//' @param undirected logical; If \code{TRUE}, only the lower half of the 
//'   adjacency matrix is used to calculate the \code{meanAdj}, halving the 
//'   calculation time.
//' @return A single numeric value.
//' @rdname meanAdj-cpp
// [[Rcpp::export]]
NumericVector MeanAdj(
  SEXP pAdjacency, IntegerVector subsetIndices, LogicalVector undirected
) {
  XPtr<BigMatrix> xpAdj(pAdjacency);
  
  // Make sure we're not indexing out of range.
  if (is_true(any(subsetIndices <= 0)) || 
      is_true(any(subsetIndices > xpAdj->ncol())) ||
      is_true(any(subsetIndices > xpAdj->nrow()))) {
    throw std::out_of_range("Requested index outside of range!");
  }
  
  // Generate appropriate warnings
  if (undirected.size() > 1) {
    Function warning("warning");
    warning("The argument 'undirected' has length > 1 and ony the first ",
            "element will be used.");
  }
  
  //  Dispatch function for all types of big.matrix.
  unsigned short type = xpAdj->matrix_type();
  if (type == 1) {
    return MeanAdj(xpAdj, MatrixAccessor<char>(*xpAdj), subsetIndices, undirected);
  } else if (type == 2) {
    return MeanAdj(xpAdj, MatrixAccessor<short>(*xpAdj), subsetIndices, undirected);
  } else if (type == 4) {
    return MeanAdj(xpAdj, MatrixAccessor<int>(*xpAdj), subsetIndices, undirected);
  } else if (type == 8) {
    return MeanAdj(xpAdj, MatrixAccessor<double>(*xpAdj), subsetIndices, undirected);
  } else {
    /* We should never get here, unless the underlying implementation of 
    bigmemory changes */
    throw Rcpp::exception("Undefined type for provided big.matrix");
  }
}
