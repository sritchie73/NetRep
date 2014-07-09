// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

/* Implementation of Intramodular Connectivity
 * 
 * @param xpAdj External Pointer for the adjacency matrix.
 * @param adj Matrix Accessor for the adjacency matrix.
 * @param subsetIndices indices of the network subset to compute the mean 
 *   adjacency on
 * @return a length one NumericVector
 */
template <typename T>
NumericVector KIM(XPtr<BigMatrix> xpAdj, MatrixAccessor<T> adj, 
                  IntegerVector subsetIndices) {
  // Results vector
  NumericVector kIM = NumericVector(subsetIndices.size());
  
  // temporary value holder
  NumericVector value = NumericVector(1); 
  
  // Make sure we're not indexing out of range.
  if (is_true(any(subsetIndices <= 0)) || 
      is_true(any(subsetIndices > xpAdj->ncol())) ||
      is_true(any(subsetIndices > xpAdj->nrow()))) {
    throw std::out_of_range("Requested index outside of range!");
  }

  int subsetSize = subsetIndices.size();

  for (int i = 0; i < subsetSize; i++) {
    kIM[i] = 0;
    for (int j = 0; j < subsetSize; j++) {
      value = adj[subsetIndices[i]-1][subsetIndices[j]-1];
      if (all(!is_na(value))) {
        kIM[i] += value[0]; // Ignore NAs
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
// [[Rcpp::export]]
NumericVector KIM(SEXP pAdjacency, IntegerVector subsetIndices) {
  //  Dispatch function for all types of big.matrix.
  XPtr<BigMatrix> xpAdj(pAdjacency);
  switch(xpAdj->matrix_type()) {
    case 1:
      return KIM(xpAdj, MatrixAccessor<char>(*xpAdj), subsetIndices);
    case 2:
      return KIM(xpAdj, MatrixAccessor<short>(*xpAdj), subsetIndices);
    case 4:
      return KIM(xpAdj, MatrixAccessor<int>(*xpAdj), subsetIndices);
    case 8:
      return KIM(xpAdj, MatrixAccessor<double>(*xpAdj), subsetIndices);
    default:
      /* We should never get here, unless the underlying implementation of 
         bigmemory changes */
      throw Rcpp::exception("Undefined type for provided big.matrix");
  }          
}
