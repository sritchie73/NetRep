// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>

/* Implementation of Mean Adjacency
 * 
 * @param xpAdj External Pointer for the adjacency matrix.
 * @param adj Matrix Accessor for the adjacency matrix.
 * @param subsetIndices indices of the network subset to compute the mean 
 *   adjacency on
 * @return a length one NumericVector
 */
template <typename T>
NumericVector MeanAdj(XPtr<BigMatrix> xpAdj, MatrixAccessor<T> adj, 
                      IntegerVector subsetIndices) {
  NumericVector mean = NumericVector(1);  // Return scalar
  
  // get some useful values
  int subsetSize = subsetIndices.size();
  int ncol = xpAdj->ncol();
  int nrow = xpAdj->nrow();
  
  // Intermediate counters
  int NAcount = 0;
  double total = 0.0;
  double value;
  
  // Make sure we're not indexing out of range.
  if (is_true(any(subsetIndices <= 0)) || is_true(any(subsetIndices > nrow)) ||
      is_true(any(subsetIndices > ncol))) {
    throw std::out_of_range("Requested index outside of range!");
  }
  
  // Add to the total sum while handling NAs
  for (int i = 0; i < subsetSize; i++) {
    for (int j = 0; j < subsetSize; j++) {
     value = adj[subsetIndices[i]-1][subsetIndices[j]-1];
      if (R_IsNA(value)) {
        NAcount += 1;
      } else { 
        total += value;
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
//' @return A single numeric value.
// [[Rcpp::export]]
NumericVector MeanAdj(SEXP pAdjacency, IntegerVector subsetIndices) {
  //  Dispatch function for all types of big.matrix.
  XPtr<BigMatrix> xpAdj(pAdjacency);
  switch(xpAdj->matrix_type()) {
    case 1:
      return MeanAdj(xpAdj, MatrixAccessor<char>(*xpAdj), subsetIndices);
    case 2:
      return MeanAdj(xpAdj, MatrixAccessor<short>(*xpAdj), subsetIndices);
    case 4:
      return MeanAdj(xpAdj, MatrixAccessor<int>(*xpAdj), subsetIndices);
    case 8:
      return MeanAdj(xpAdj, MatrixAccessor<double>(*xpAdj), subsetIndices);
    default:
      /* We should never get here, unless the underlying implementation of 
         bigmemory changes */
      throw Rcpp::exception("Undefined type for provided big.matrix");
  }   
}
