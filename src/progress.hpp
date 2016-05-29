#ifndef __PROGRESS__
#define __PROGRESS__

#include <RcppArmadillo.h>
#include "interrupt.hpp"

void MonitorProgress (unsigned int&, arma::uvec&, bool&, const bool&); 

#endif // __PROGRESS__
