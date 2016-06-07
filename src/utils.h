#ifndef __UTILS__
#define __UTILS__

#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#define ARMA_NO_DEBUG
#define ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_USE_CXX11
#define BOOST_DISABLE_ASSERTS

#include <RcppArmadillo.h>
#include <boost/unordered_map.hpp>
#include <string>

// For mapping column/row/names to indices of respective data structures
typedef boost::unordered_map<std::string, unsigned int> namemap; 
// For mapping module labels to node IDs 
typedef boost::unordered_multimap<std::string, std::string> stringmap; 
// For storing property vectors across modules
typedef boost::unordered_map<std::string, arma::vec> vecmap;
// For storing their addresses to pass to threads
typedef boost::unordered_map<std::string, double *> addrmap;
// For getting randomly shuffled node ids
typedef boost::unordered_map<unsigned int, unsigned int> intmap;

// Utility functions
void ShowProgress(unsigned int&, unsigned int&);
namemap MakeIdxMap (const std::vector<std::string>&);
stringmap MakeModMap (Rcpp::CharacterVector);
stringmap MakeModMap (Rcpp::CharacterVector, const namemap&);
namemap MakeNullMap (const std::vector<std::string>&, const namemap&, arma::uvec&);
arma::uvec GetNodeIdx (std::string&, const stringmap&, const namemap&);
arma::uvec GetRandomIdx(std::string&, const stringmap&, arma::uvec&, namemap&);
std::vector<std::string> GetModNodeNames (std::string&, const stringmap&);
void Fill(Rcpp::NumericVector&, arma::vec&, arma::uvec&);

#endif // __UTILS__