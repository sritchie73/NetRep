#include "utils.h"

/* Build a dictionary mapping labels to a sequence of integers
 */
namemap MakeIdxMap (const std::vector<std::string>& labels) {
  namemap map;
  for (unsigned int ii = 0; ii < labels.size(); ++ii) {
    map[labels[ii]] = ii;
  }
  return map;
}

/* Builds the dictionary between module labels and their associated nodes
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

/* Builds the dictionary between module labels and their associated nodes,
 * containing only nodes that are present in the dataset of interest.
 *
 * @param moduleAssignments a named vector of module labels where the names are
 *   the node IDs.
 * @param nodeIdx a mapping of node labels to indices in the dataset of interest.
 *
 * @return an unordered_map where the keys are module labels, and values are
 *   arrays containing the node IDs belonging to that module
 */
stringmap MakeModMap (
    Rcpp::CharacterVector moduleAssignments, const namemap& nodeIdx
) {
  stringmap map;
  
  const std::vector<std::string> nodes = moduleAssignments.names();
  const std::vector<std::string> labels = Rcpp::as<std::vector<std::string>>(moduleAssignments);
  
  for (unsigned int ii = 0; ii < moduleAssignments.length(); ++ii) {
    std::string key = labels[ii];
    std::string value = nodes[ii];
    
    if (nodeIdx.find(value) != nodeIdx.end()) {
      map.emplace(key, value);
    }
  }
  
  return map;
}

/* Build the vector and map for use when generating null distributions
 * 
 * The permutation procedure has two options: (1) 'nullHypothesis' == "all",
 * in which case the null distributions are generated using all nodes in the
 * test network as a background, or (2) 'nullHypothesis' == "overlap", in
 * which case the null distributions are generated using only nodes present
 * in both the discovery and test datasets (this is appropriate where the 
 * assumption is that unobserved nodes from the test dataset but not present
 * in the discovery dataset may very well be present in the module of 
 * interest).
 * 
 * Internally, the C++ code works by iterating through each module of interest,
 * then pulling out the node labels associated with that module from a map 
 * (O(1) lookup) then pulling out the indices in each dataset associated with
 * those nodes from another map (O(1) lookup). Since random number generation 
 * is non-trivial, and computationally expensive, the simplest solution is to
 * shuffle the test dataset at each permutation. To do this, we will have a 
 * vector of node IDs that can be used when generating the null hypothesis, and
 * a mapping between node IDs and indices in that vector. Then, at each 
 * permutation, we can shuffle the 'nullIdx', then the corresponding lookups 
 * will go module -> node IDs (O(1) lookup), node IDs -> index in 'nullMap' 
 * (O(1) lookup), and then index in 'nullMap' -> random indice in the test 
 * dataset.
 * 
 * @param validNodes a vector of node names to be used when generating the 
 *   null hypothesis.
 * @param tIdxMap a mapping of node IDs to indices in the test dataset.
 * @param nullIdx a vector that has been initialised to the size of 
 *   'validNodes', this will be filled with indices in the test dataset.
 * 
 * @details
 *   If the node names of the test dataset are provided to 'validNodes', then
 *   all nodes in the test dataset will be used (i.e. 'nullHypothesis' = "all").
 *   If the node names of the discovery dataset are provided, then only nodes
 *   that are present in both the discovery and test dataset will be used (i.e
 *   'nullHypothesis' = "overlap").
 * 
 * @return
 *   a mapping from node IDs to their indice in the filled 'nullIdx'. 
 *   'nullIdx' will be filled with valid indices as a side effect.
 */
namemap MakeNullMap (
  const std::vector<std::string>& validNodes, const namemap& tIdxMap, 
  arma::uvec& nullIdx
) {
  namemap nullMap;
  nullIdx.set_size(validNodes.size()); // Will have, at most, this many elements
  
  unsigned int counter = 0;
  // For each valid node
  for (unsigned int ii = 0; ii < validNodes.size(); ++ii) {
    // If the node is present in the test dataset
    if (tIdxMap.count(validNodes[ii]) == 1) {
      // Insert its position in the test dataset into our vector of valid null
      // indices
      nullIdx.at(counter) = tIdxMap.at(validNodes[ii]);
      // And tell our 'nullMap' about the position of this node in 'nullIdx'
      nullMap[validNodes[ii]] = counter;
      counter++;
    }
  }
  
  // Now we can shrink 'nullIdx' to the total number of null nodes in the case
  // of 'nullHypothesis' == "overlap"
  if (counter < validNodes.size()) {
    nullIdx.resize(counter);
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

/* Get the indices of a module's nodes in the respective dataset
 * 
 * @param module module we want to get the indices for
 * @param modNodeMap mapping between module labels and node IDs.
 * @param nodeIdxMap mapping between node IDs and indices in the dataset of
 *   interest
 *
 * @return a vector of indices   
 */
std::vector<std::string> GetModNodeNames (
    std::string& module, const stringmap& modNodeMap
) {
  unsigned int nNodes = modNodeMap.count(module); 
  std::vector<std::string> modNames (nNodes);
  
  unsigned int counter = 0;
  auto keyit = modNodeMap.equal_range(module);
  for (auto it = keyit.first; it != keyit.second; ++it) {
    std::string nodeID = it->second;
    modNames.at(counter) = nodeID;
    counter++;
  }
  
  return modNames;
}

/* Fill a non-contiguous subset of a NumericVector with the contents of an
 * armadillo vector
 * 
 * Note: 'idx' must be the same size as 'contents': this is not checked by the 
 * code!
 *
 * @param tofill vector to fill.
 * @param contents vector whose contents to fill 'tofill' with.
 * @param idx (ordered) indices in 'tofill' to use when filling 'tofill'.
 * 
 */
void Fill(Rcpp::NumericVector& tofill, arma::vec& contents, arma::uvec& idx) {
  for (unsigned int ii=0; ii < contents.size(); ++ii) {
    tofill[idx.at(ii)] = contents.at(ii);
  }
}
