/**
 * This file contains utility functions for dealing with big.matrix objects
 */

// For safely extracting single values from a big.matrix.
NumericVector safeAccessor(XPtr<BigMatrix> pBigMat, int row, int column);
