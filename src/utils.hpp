#ifndef __UTILS__
#define __UTILS__

#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>

// For mapping column/row/names to indices of respective data structures
typedef std::unordered_map<std::string, unsigned int> namemap; 
// For mapping module labels to node IDs 
typedef std::unordered_multimap<std::string, std::string> stringmap; 
// For storing properties across modules
typedef std::unordered_map<std::string, arma::vec> vecmap;

// Utility functions
namemap makeIdxMap (std::vector<std::string>);
stringmap MakeModMap (Rcpp::CharacterVector);
arma::uvec GetNodeIdx (std::string&, const stringmap&, const namemap&);

#endif // __UTILS__