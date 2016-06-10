#include "scale.h"

/* Scale data across all nodes
 * 
 * Each node is centered by its mean and scaled by it standard deviation.
 * 
 * @param dataAddr memory address of the scaled data matrix.
 * @param nSamples number of samples in the dataset.
 * @param nNodes number of nodes in the network.
 * 
 * @return
 *  A scaled data matrix.
 */
arma::mat Scale (double * dataAddr, unsigned int nSamples, unsigned int nNodes) {
  // Construct armadillo matrices and vectors from memory addresses and 
  // provided sizes
  arma::mat data = arma::mat(dataAddr, nSamples, nNodes, false, true);
  arma::mat scaled = arma::mat(data.n_rows, data.n_cols);
  
  for (unsigned int ii = 0; ii < data.n_cols; ++ii) {
    arma::vec nodeData = data.col(ii);
    scaled.col(ii) = (nodeData - mean(nodeData))/stddev(nodeData, 0);
  }
  
  return scaled;
}

//' Scale data across all nodes
//' 
//' Each node is centered by its mean and scaled by it standard deviation.
//' 
//' @param data matrix to scale.
//' 
//' @return
//'  A scaled data matrix.
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix Scale(Rcpp::NumericMatrix data) {
  arma::mat scaled = Scale(data.begin(), data.nrow(), data.ncol());
  Rcpp::NumericMatrix forR (scaled.n_rows, scaled.n_cols, scaled.begin());
  colnames(forR) = colnames(data);
  rownames(forR) = rownames(data);
  return(forR);
}
