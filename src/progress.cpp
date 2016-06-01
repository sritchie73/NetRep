#if WINDOWS
  #include <windows.h>
#else
  #include <unistd.h>
#endif

#include "progress.hpp"


// Sleep this thread for N seconds
void sleep(int secs) {
  #if WINDOWS
    Sleep(secs * 1000);
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
 * @param progress a vector that is being updated by each thread with the 
 *   number of permutations that have been completed by each thread
 * @param interrupted a boolean value in the heap, accessible to each thread,
 *   that is modified when the user attempts to cancel the permutation 
 *   procedure in the R session.
 * @param verboseFlag if 'false' messages are not printed.
 */
void MonitorProgress (
    unsigned int& nPerm, arma::uvec& progress, bool& interrupted,
    const bool& verboseFlag
) {
  if (verboseFlag) {
    Rcpp::Rcout << "\n";
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
      sprintf(formatted, "%5d", percentCompleted);
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
    Rcpp::Rcout << "\n\n";
  }
}
