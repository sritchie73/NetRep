#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#define ARMA_NO_DEBUG
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>

// header files for sleeping a thread
#if defined (_WIN32) || defined (_WIN64)
  // <thread> requires Rtools33 (R 3.3.0) which uses GCC 4.9.3, which supports 
  // <chrono>
  #include <chrono>
#else
  // But for linux/osx we need to handle older versions of GCC which have 
  // <thread> but not <chrono> support
  #include <unistd.h>
#endif

#include "thread-utils.h"

// Sleep this thread for N seconds
void sleep(int secs) {
  #if defined (_WIN32) || defined (_WIN64)
    std::this_thread::sleep_for(std::chrono::seconds(secs));
  #else
    usleep(secs * 1000000);  
  #endif
}

/* Monitors the progress of the permutation procedure
 * 
 * The primary role of this function is to facilitate user interrupts to cancel
 * the permutation procedure. This runs every second until the permutaton 
 * procedure is complete. If 'verboseFlag' is 'true' then the current progress 
 * is also printed every second.
 * 
 * @param nPerm the total number of permutations to compute
 * @param progressAddr memory address of a vector that is being updated by each 
 *   thread with the number of permutations that have been completed by each 
 *   thread.
 * @param nThreads total number of threads being run.
 * @param interrupted a boolean value in the heap, accessible to each thread,
 *   that is modified when the user attempts to cancel the permutation 
 *   procedure in the R session.
 * @param verboseFlag if 'false' messages are not printed.
 */
void MonitorProgress (
    unsigned int& nPerm, unsigned int * progressAddr, unsigned int nThreads,
    bool& interrupted, const bool& verboseFlag
) {
  arma::uvec progress = arma::uvec(progressAddr, nThreads, false, true);
  
  if (verboseFlag) {
    Rcpp::Rcout << std::endl;
  }
  unsigned int nCompleted = 0;
  unsigned int percentCompleted = 0;
  char formatted[6]; // stores a whitespace padded percentage value
  
  // Every second, check the number of permutations completed, print out the 
  // % completed, and check for interrupts from the R session.
  while (true) {
    nCompleted = sum(progress);
    if (verboseFlag) {
      percentCompleted = (unsigned int)round( (float)nCompleted / (float)nPerm * 100);
      snprintf(formatted, sizeof(formatted), "%5d", percentCompleted);
      Rcpp::Rcout << "\r" << formatted << "% completed."; 
    }
    if (nCompleted == nPerm) {
      break;
    } 
    if (checkInterrupt()) {
      interrupted = true;
      break;
    }
    sleep(1);
  }
  if (verboseFlag) {
    Rcpp::Rcout << std::endl << std::endl;
  }
}
