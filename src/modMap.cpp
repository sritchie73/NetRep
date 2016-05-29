#include "utils.hpp"

#define ARMA_NO_DEBUG

/* Build a dictionary mapping labels to a sequence of integers
 */
namemap MakeIdxMap (const std::vector<std::string>& labels) {
  namemap map;
  for (unsigned int ii = 0; ii < labels.size(); ++ii) {
    map[labels[ii]] = ii;
  }
  return map;
}

/* Builds the dictionary between module labels are their associated nodes
 * 
 * @param moduleAssignments a named vector of module labels where the names are
 *   the node IDs.
 *
 * @return an unordered_map where the keys are module labels, and values are
 *   arrays containing the node IDs belonging to that module
 */
stringmap MakeModMap (Rcpp::CharacterVector moduleAssignments) {
  stringmap map;
  
  const std::vector<std::string> nodes = moduleAssignments.names();
  const std::vector<std::string> labels = Rcpp::as<std::vector<std::string>>(moduleAssignments);
  
  for (unsigned int ii = 0; ii < moduleAssignments.length(); ++ii) {
    std::string key = labels[ii];
    std::string value = nodes[ii];
    
    map.emplace(key, value);
  }
  
  return map;
}

/* Build a vector of in 
 * 
 * @param validNodes a vector of node IDs to use when generating the null 
 *   distributions.
 * @param tIdxMap a mapping of node IDs to indices in the test dataset.
 * @param nullIdx a vector that has been initialised to the size of 
 *   'validNodes', this will be filled with indices in the test dataset.
 * 
 * @return
 *   a mapping from node IDs to their indice in the filled 'nullIdx'. 
 *   'nullIdx' will be filled with valid indices as a side effect.
 */
namemap MakeNullMap (
  const std::vector<std::string>& validNodes, const namemap& tIdxMap, 
  arma::uvec& nullIdx
) {
   (validNodes.size());
  namemap nullMap;
  
  // For each valid node
  for (unsigned int ii = 0; ii < validNodes.size(); ++ii) {
    // Look up its index in the test dataset, and insert it into the nullNodes
    // vector, which will have the same ordering of nodes as validNodes.
    nullIdx.at(ii) = tIdxMap.at(validNodes[ii]);
    // Insert its index in 'nullNodes' into the map 
    nullMap[validNodes[ii]] = ii;
  }
  
  return nullMap;
}

/* Get the indices of a module's nodes in the respective dataset
 * 
 * @param module module we want to get the indices for
 * @param modNodeMap mapping between module labels and node IDs.
 * @param nodeIdxMap mapping between node IDs and indices in the dataset of
 *   interest
 *
 * @return a vector of indices   
 */
arma::uvec GetNodeIdx (
  std::string& module, const stringmap& modNodeMap, const namemap& nodeIdxMap
) {
  unsigned int nNodes = modNodeMap.count(module); 
  arma::uvec modIdx (nNodes);
  
  unsigned int counter = 0;
  auto keyit = modNodeMap.equal_range(module);
  for (auto it = keyit.first; it != keyit.second; ++it) {
    std::string nodeID = it->second;
    modIdx.at(counter) = nodeIdxMap.at(nodeID);
    counter++;
  }

  return modIdx;
}

/* Get a random selection of nodes of size N from a dataset
 * 
 * @param module module we want to get random indices for.
 * @param modNodeMap mapping between module labels and node IDs.
 * @param nodeIdx a (shuffled) vector of node indices in the test dataset. 
 *  These indices correspond to the set of nodes to use when generating the 
 *  null distributions.
 * @param nullMap a mapping of node IDs to indices of the 'nodeIdx' vector.
 * 
 * @return a vector of indices in the test dataset.
 */ 
arma::uvec GetRandomIdx(
  std::string& module, const stringmap& modNodeMap, arma::uvec& nodeIdx, 
  namemap& nullMap
) {
  unsigned int nNodes = modNodeMap.count(module);
  arma::uvec randIdx (nNodes);
  
  // For each node in the module, get the nodes static position in the 'nodeIdx'
  // vector, pull out the randomly assigned indice in the test network 
  // stored in that location in 'nodeIdx', and put it in the vector of random
  // node ids. Random assignment of indices happens once per permutation, 
  // through the use of 'shuffle' on the 'nodeIdx' vector.
  unsigned int counter = 0;
  auto keyit = modNodeMap.equal_range(module);
  for (auto it = keyit.first; it != keyit.second; ++it) {
    std::string nodeId = it->second;
    randIdx.at(counter) = nodeIdx.at(nullMap.at(nodeId));
    counter++;
  }
  
  return randIdx;
}
