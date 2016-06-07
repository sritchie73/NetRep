#ifndef __PROGRESS__
#define __PROGRESS__

#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#define ARMA_NO_DEBUG
#define ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_USE_CXX11

#include <RcppArmadillo.h>
#include "interrupt.h"
#include <thread>

void MonitorProgress (unsigned int&, arma::uvec&, bool&, const bool&); 

#endif // __PROGRESS__
