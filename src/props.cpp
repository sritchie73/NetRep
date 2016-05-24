/* Functions for calculating the network properties (but not statistics) for a
 * single module within a single dataset.
 */
#include "functions.hpp"

// Sort nodes for sequential memory access
//
// When running calculations on an arbitrary sub-matrix (e.g. a module, or a
// random sub-graph during the permutation procedure) it is much more efficient
// to sort the requested indices so that memory access is sequential wherever
// possible. However, we need to re-order the results vector so that the
// properties for each node are in the same order as requested in the function.
//
// @return
//   Sorts 'nodeIdx' as a side effect and returns a vector of ranks that can
//   be used to re-order another results vector so that nodes are in the same
//   order as in 'nodeIdx' prior to sorting.
uvec sortNodes (uvec& nodeIdx) {
  uvec rank = arma::regspace<uvec>(0, nodeIdx.n_elem);  // seq_along(nodeIdx)
  uvec order = sort_index(nodeIdx);
  nodeIdx = nodeIdx(order);
  return rank(order);
}

// Calculate the weighted degree of a module
//
// The weighted degree is the sum of edge weights to all other nodes in the
// network. Assumes that the network is undirected.
//
// @param netPtr address of the network's adjacency matrix in memory.
// @param nodeIdx ascending-order sorted indices of the module's nodes.
//
// @return a (column) vector of the weighted degree.
//
vec WeightedDegree(const mat& netPtr, uvec& nodeIdx) {
  // We do not want a negative weight to cancel out a positive one, so we take
  // the absolute value.
  mat dg = diagvec(netPtr);
  mat wDegree = sum(abs(netPtr(nodeIdx, nodeIdx))) - abs(dg(nodeIdx)).t();

  return wDegree;
}

// Calculate the average edge weight
//
// @param wDegree the weighted degree, see 'WeightedDegree'
//
// @return a scalar value
double AverageEdgeWeight(vec& wDegree) {
  int n = wDegree.n_elem;
  return as_scalar(sum(wDegree) / (n*n - n));
}

// Get a vector of correlation coefficients for a module
//
// Flattens a sub-matrix view of the corrPtr matrix minus the diagonals
//
// @param corrPtr address in memory of the matrix of correlation coefficients.
// @param nodeIdx ascending-order sorted indices of the module's nodes.
//
// @return a vector of correlation coefficients
vec CorrVector (const mat& corrPtr, uvec& nodeIdx) {
  // Number of nodes in the requested sub-matrix
  int n = nodeIdx.n_elem;
  
  // We need to flatten the matrices to a vector, ignoring the diagonals.
  unsigned int flatsize = (n*n - n)/2;
  vec corrVec(flatsize);
  
  unsigned int vi = 0; // keeps track of position in corrVec
  
  // Iterate over columns and rows to fill out 'corrVec' with the lower 
  // triangle of the submatrix.
  for (unsigned int jj = 0; jj < n; jj++) {
    for (unsigned int ii = jj + 1; ii < n; ii++) {
      corrVec.at(vi) = corrPtr(nodeIdx(ii), nodeIdx(jj));
      vi++;
    }   
  }
  
  return corrVec;
}

// Calculate the summary profile of a module
// 
// @param dataPtr address of the data matrix in memory
// @param nodeIdx ascending-order sorted indices of the module's nodes.
// 
// @return a vector of observations across samples
vec SummaryProfile (const mat& dataPtr, uvec& nodeIdx) {
  mat U, V;
  vec S;
  
  bool success = svd_econ(U, S, V, dataPtr.cols(nodeIdx), "left", "dc");
  
  if (!success) {
    vec summary (1);
    summary.fill(-9);
  }
  
  vec summary = U.col(0);
  
  // Flip the sign of the summary profile so that the eigenvector is 
  // positively correlated with the average scaled value of the underlying
  // data for the network module.
  vec meanObs = mean(dataPtr.cols(nodeIdx), 1);
  int orientation = as_scalar(sign(cor(meanObs, summary)));
  
  if (orientation == -1) {
    summary *= -1; 
  }
  
  return summary;
}

// Calculate the contribution of each node to the summary profile
//
// @param dataPtr address in memory of the data matrix.
// @param nodeIdx ascending-order sorted indices of the module's nodes.
// @param summaryProfile the summary profile vector, see 'SummaryProfile'
//
// @return a vector of correlations between each node and the summary profile
vec NodeContribution (const mat& dataPtr, uvec& nodeIdx, vec& summaryProfile) {
  return cor(summaryProfile, dataPtr.cols(nodeIdx));
}

// Calculate module's coherence
// 
// As measured by the proportion of variance in the module's data explained 
// by the module summary profile.
//
// @param nodeContribution the vector of node contributions, see 'NodeContribution'
//
// @return a double between 0 and 1
double ModuleCoherence (vec& nodeContribution) {
  return mean(square(nodeContribution));
}


