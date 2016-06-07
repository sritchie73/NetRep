/* Functions for calculating the network properties and module preservation 
 * statistics for a single module within a single dataset.
 */

#include "netStats.h"

/* Scale data across all nodes
 * 
 * Each node is centered by its mean and scaled by it standard deviation.
 */
arma::mat Scale (const arma::mat& dataPtr) {
  arma::mat scaled = arma::mat(dataPtr.n_rows, dataPtr.n_cols);
  
  for (unsigned int ii = 0; ii < dataPtr.n_cols; ++ii) {
    arma::vec nodeData = dataPtr.col(ii);
    scaled.col(ii) = (nodeData - mean(nodeData))/stddev(nodeData, 0);
  }
  
  return scaled;
}

/* Sort nodes for sequential memory access
 *
 * When running calculations on an arbitrary sub-matrix (e.g. a module, or a
 * random sub-graph during the permutation procedure) it is much more efficient
 * to sort the requested indices so that memory access is sequential wherever
 * possible. However, we need to re-order the results vector so that the
 * properties for each node are in the same order as requested in the function.
 *
 * @return
 *   Sorts 'nodeIdx' as a side effect and returns a vector of ranks that can
 *   be used to re-order another results vector so that nodes are in the same
 *   order as in 'nodeIdx' prior to sorting.
 */
arma::uvec sortNodes (arma::uvec& nodeIdx) {
  arma::uvec rank = arma::sort_index(arma::sort_index(nodeIdx));
  arma::uvec order = arma::sort_index(nodeIdx);
  nodeIdx = nodeIdx(order);
  return rank;
}

/* Calculate the correlation between two vectors
* @param v1addr memory address of a vector.
* @param v2addr memory address of a vector.
* @param size size of the two vectors
*/
double Correlation (double * v1addr, double * v2addr, unsigned int size) {
  arma::vec v1 = arma::vec(v1addr, size, false, true);
  arma::vec v2 = arma::vec(v2addr, size, false, true);
  
  return arma::as_scalar(arma::cor(v1, v2));
}

/* Calculate the sign-aware mean of two vectors
 * 
 * This is the mean of 'v2' where observations detract from the mean if they
 * differ in sign between 'v1' and 'v2'
 * 
 * @param v1addr memory address of a vector.
 * @param v2addr memory address of a vector.
 * @param size size of the two vectors
 * 
 */
double SignAwareMean (double * v1addr, double * v2addr, unsigned int size) {
  arma::vec v1 = arma::vec(v1addr, size, false, true);
  arma::vec v2 = arma::vec(v2addr, size, false, true);
  return arma::mean(arma::sign(v1) % v2);
}

/* Calculate the weighted degree of a module
 *
 * The weighted degree is the sum of edge weights to all other nodes in the
 * network. Assumes that the network is undirected.
 *
 * @param netPtr address of the network's adjacency matrix in memory.
 * @param nodeIdx ascending-order sorted indices of the module's nodes.
 *
 * @return a (column) vector of the weighted degree.
 */
arma::vec WeightedDegree(const arma::mat& netPtr, arma::uvec& nodeIdx) {
  // We take the absolute value so that negative weights (if they exist) do not
  // cancel out positive ones
  arma::rowvec colSums = arma::sum(arma::abs(netPtr(nodeIdx, nodeIdx)), 0);
  // We need to convert to a column-vector
  // Note to self: do not set copy_aux_memory to false! The program may free
  // the memory used in the returned vector to be used elsewhere!
  arma::vec wDegree = arma::vec(colSums.begin(), colSums.n_elem, true);
  // subtract the diagonals
  arma::vec dg = netPtr.diag(0);
  wDegree -= arma::abs(dg.elem(nodeIdx));
  return wDegree;
}

/* Calculate the average edge weight
 *
 * @param wDegreeAddr memory address of the weighted degree vector, see 
 *  'WeightedDegree'.
 * @param nNodes number of nodes in the module
 *
 * @return a scalar value
 */
double AverageEdgeWeight(double * wDegreeAddr, unsigned int nNodes) {
  arma::vec wDegree = arma::vec(wDegreeAddr, nNodes, false, true);
  double nEdgePairs = (double)(nNodes*nNodes - nNodes);
  double allEdges =  arma::as_scalar(arma::sum(wDegree));
  return allEdges / nEdgePairs;
}

/* Get a vector of correlation coefficients for a module
 *
 * Flattens a sub-matrix view of the corrPtr matrix minus the diagonals
 *
 * @param corrPtr address in memory of the matrix of correlation coefficients.
 * @param nodeIdx ascending-order sorted indices of the module's nodes.
 *
 * @return a vector of correlation coefficients
 */
arma::vec CorrVector (const arma::mat& corrPtr, arma::uvec& nodeIdx) {
  // Number of nodes in the requested sub-matrix
  int n = nodeIdx.n_elem;
  
  // We need to flatten the matrices to a vector, ignoring the diagonals.
  unsigned int flatsize = (n*n - n)/2;
  arma::vec corrVec(flatsize);
  
  unsigned int vi = 0;  // keeps track of position in corrVec
  
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

/* Calculate the summary profile of a module
 * 
 * @param dataPtr address of the data matrix in memory
 * @param nodeIdx ascending-order sorted indices of the module's nodes.
 * 
 * @return a vector of observations across samples
 */
arma::vec SummaryProfile (const arma::mat& dataPtr, arma::uvec& nodeIdx) {
  arma::mat U, V;
  arma::vec S;
  
  bool success = arma::svd_econ(U, S, V, dataPtr.cols(nodeIdx), "left", "dc");
  
  if (!success) {
    arma::vec summary (1);
    summary.fill(arma::datum::nan);
    return summary;
  }
  arma::vec summary = U.col(0);
  
  /* Flip the sign of the summary profile so that the eigenvector is 
   * positively correlated with the average scaled value of the underlying
   * data for the network module.
   */
  arma::vec meanObs = arma::mean(dataPtr.cols(nodeIdx), 1);
  int orientation = arma::as_scalar(arma::sign(arma::cor(meanObs, summary)));

  if (orientation == -1) {
    summary *= -1;
  }
  
  return summary;
}

/* Calculate the contribution of each node to the summary profile
 *
 * @param dataPtr address in memory of the data matrix.
 * @param nodeIdx ascending-order sorted indices of the module's nodes.
 * @param summaryProfile the summary profile vector, see 'SummaryProfile'
 *
 * @return a vector of correlations between each node and the summary profile
 */
arma::vec NodeContribution (
  const arma::mat& dataPtr, arma::uvec& nodeIdx, arma::vec& summaryProfile
) {
  // We need to convert SP from a vector to a matrix since arma::cor doesn't
  // have a method for comparing a matrix to a vector.
  const arma::mat SP = arma::mat(summaryProfile.begin(), summaryProfile.n_elem, 1, true);
  return arma::cor(dataPtr.cols(nodeIdx), SP);
}

/* Calculate module's coherence
 * 
 * As measured by the proportion of variance in the module's data explained 
 * by the module summary profile.
 *
 * @param ncAddr memory address of the vector of node contributions, see 
 *  'NodeContribution'.
 * @param nNodes number of nodes in the module.
 *
 * @return a double between 0 and 1
 */
double ModuleCoherence (double * ncAddr, unsigned int nNodes) {
  arma::vec nodeContribution = arma::vec(ncAddr, nNodes, false, true);
  return arma::mean(arma::square(nodeContribution));
}


