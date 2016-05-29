// Functions for checking whether the user has asked the C++ code to interrupt

#include <RcppArmadillo.h>

using namespace Rcpp;

static void chkIntFn(void *dummy) { 
  R_CheckUserInterrupt(); 
} 

// this will call the above in a top-level context so it won't longjmp-out of your context 
bool checkInterrupt() { 
  return (R_ToplevelExec(chkIntFn, NULL) == FALSE); 
} 
