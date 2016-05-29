#ifndef __UTILS__
#define __UTILS__

#include <RcppArmadillo.h>
#include <boost/unordered_map.hpp>
#include <string>

// For mapping column/row/names to indices of respective data structures
typedef boost::unordered_map<std::string, unsigned int> namemap; 
// For mapping module labels to node IDs 
typedef boost::unordered_multimap<std::string, std::string> stringmap; 
// For storing properties across modules
typedef boost::unordered_map<std::string, arma::vec> vecmap;
// For getting randomly shuffled node ids
typedef boost::unordered_map<unsigned int, unsigned int> intmap;

// Utility functions
void ShowProgress(unsigned int&, unsigned int&);
namemap MakeIdxMap (const std::vector<std::string>&);
stringmap MakeModMap (Rcpp::CharacterVector);
namemap MakeNullMap (const std::vector<std::string>&, const namemap&, arma::uvec&);
arma::uvec GetNodeIdx (std::string&, const stringmap&, const namemap&);
arma::uvec GetRandomIdx(std::string&, const stringmap&, arma::uvec&, namemap&);

#endif // __UTILS__