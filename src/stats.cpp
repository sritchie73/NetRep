/* Functions for calculating the network properties (but not statistics) for a
 * single module within a single dataset.
 */
#include "functions.hpp"

// Calculate the correlation between two vectors
double Correlation (vec& v1, vec& v2) {
  return as_scalar(cor(v1, v2));
}

// Calculate the sign-aware mean of two vectors
// 
// This is the mean of 'v2' where observations detract from the mean if they
// differ in sign between 'v1' and 'v2'
//
double SignAwareMean (vec& v1, vec& v2) {
  return mean(sign(v1) % v2);
}