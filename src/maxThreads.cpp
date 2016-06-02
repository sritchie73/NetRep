#include <Rcpp.h>
#include <thread>

//' Get the maximum number of concurrent threads supported by this machine
//' 
// [[Rcpp::export]]
Rcpp::NumericVector MaxThreads () {
  return Rcpp::wrap( std::thread::hardware_concurrency() );
}
