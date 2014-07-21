#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

#include <cmath>

/* Implementation of Intramodular Connectivity
 * 
 * @param xpAdj External Pointer for the adjacency matrix.
 * @param adj Matrix Accessor for the adjacency matrix.
 * @param subsetIndices indices of the network subset to compute the mean 
 *   adjacency on
 * @return a length one NumericVector
 */
template <typename T>
NumericVector KIM(
  XPtr<BigMatrix> xpAdj, MatrixAccessor<T> adj, IntegerVector subsetIndices
) {
  // Results vector
  NumericVector kIM(subsetIndices.size(), 0.0);
  
  // get some useful values
  int subsetSize = subsetIndices.size();
  int type = xpAdj->matrix_type();
  
  // temporary value holder
  double value; 

  for (unsigned int jj = 0; jj < subsetSize; jj++) {
    for (unsigned int ii = 0; ii < subsetSize; ii++) {
      value = adj[subsetIndices[jj]-1][subsetIndices[ii]-1];
      if (type == 1 && !(value == NA_CHAR)) {
        kIM[jj] += abs(value);
      } else if (type == 2 && !(value == NA_SHORT)) {
        kIM[jj] += abs(value);
      } else if (type == 4 && !(value == NA_INTEGER)) {
        kIM[jj] += abs(value);
      } else if (type == 8 && !ISNA(value) && !ISNAN(value)) {
        kIM[jj] += abs(value);
      }
    }
  }
  
  return kIM;
}

//' C++ implementation of Intramodular Connectivity (Degree).
//' 
//' See the \link[=kIM]{wrapper function} for documentation.
//' 
//' @param pAdjacency SEXP container for the pointer to the adjacency matrix
//' @param subsetIndices indices of the subset of the network to calculate
//'   the mean adjacency for.
//' @return A vector containing the intramodular connectivity (degree) of 
//'   each node. 
//' @rdname kIM-cpp
// [[Rcpp::export]]
NumericVector KIM(SEXP pAdjacency, IntegerVector subsetIndices) {
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
    return KIM(xpAdj, MatrixAccessor<char>(*xpAdj), subsetIndices);
  } else if (type == 2) {
    return KIM(xpAdj, MatrixAccessor<short>(*xpAdj), subsetIndices);
  } else if (type == 4) {
    return KIM(xpAdj, MatrixAccessor<int>(*xpAdj), subsetIndices);
  } else if (type == 8) {
    return KIM(xpAdj, MatrixAccessor<double>(*xpAdj), subsetIndices);
  } else {
    /* We should never get here, unless the underlying implementation of 
    bigmemory changes */
    throw Rcpp::exception("Undefined type for provided big.matrix");
  }
}
