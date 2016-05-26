#include "utils.hpp"

using namespace std;
using namespace arma;
using namespace Rcpp;


/* Build a dictionary mapping labels to a sequence of integers
 */
namemap makeIdxMap (vector<string> labels) {
  namemap map;
  for (unsigned int ii = 0; ii < labels.size(); ++ii) {
    map.emplace(labels[ii], ii);
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
stringmap MakeModMap (CharacterVector moduleAssignments) {
  stringmap map;
  
  const vector<string> nodes = moduleAssignments.names();
  const vector<string> labels = as<vector<string>>(moduleAssignments);
  
  for (unsigned int ii = 0; ii < moduleAssignments.length(); ++ii) {
    string key = labels[ii];
    string value = nodes[ii];
    
    map.emplace(key, value);
  }
  
  return map;
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
uvec GetNodeIdx (
  string& module, const stringmap& modNodeMap, const namemap& nodeIdxMap
) {
  unsigned int nNodes = modNodeMap.count(module); 
  uvec modIdx (nNodes);
  
  unsigned int counter = 0;
  auto keyit = modNodeMap.equal_range(module);
  for (auto it = keyit.first; it != keyit.second; ++it) {
    string nodeID = it->second;
    modIdx(counter) = nodeIdxMap.at(nodeID);
    counter++;
  }

  return modIdx;
}
