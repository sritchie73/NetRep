#define ARMA_USE_LAPACK
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
#include <bigmemory/MatrixAccessor.hpp>

// Implementation of CoexpStats
//
// @param xpCoexpD,xpCoexpT external pointers to the 'big.matrix' objects for 
//  the discovery and test correlation networks respectively.
// @param matCoexpD,matCoexpT Matrix Accessors for the 'big.matrix' objects for
//  the discovery and test correlation networks respectively.
// @param dIdx,tIdx indices of the network subset in the discovery and test
//  networks respectively.
//
// @return A list containing the
//    - correlation of the correlation (cor.coexp)
//    - sign aware mean of the correlation (mean.coexp)
template <typename S, typename T>
List CoexpStats(
  XPtr<BigMatrix> xpCoexpD, MatrixAccessor<S> matCoexpD, IntegerVector dIdx,
  XPtr<BigMatrix> xpCoexpT, MatrixAccessor<T> matCoexpT, IntegerVector tIdx
) {
  // Number of nodes in the network subset
  int n = dIdx.size();

  // We need to flatten the matrices to a vector, ignoring the diagonals.
  unsigned int flatsize = (n*n - n)/2;
  rowvec corD(flatsize);
  rowvec corT(flatsize);
  unsigned int vi = 0; // keeps track of position in corD and corT.

  for (unsigned int jj = 0; jj < n; jj++) {
    for (unsigned int ii = 0; ii < jj; ii++) {
      corD.at(vi) = (double)matCoexpD[dIdx[jj] - 1][dIdx[ii] - 1]; 
      corT.at(vi) = (double)matCoexpT[tIdx[jj] - 1][tIdx[ii] - 1]; 
      vi++;
    }   
  }

  mat corCor = cor(corD, corT);
  mat meanCor = mean(sign(corD) % corT, 1); 


  if (corCor(0,0) < -1 || corCor(0,0) > 1) {
    Function warning("warning");
    warning("'cor.coexp' returned a correlation outside of [-1,1]"
            "returning NA instead.");
    return List::create(
      Named("cor.coexp") = NA_REAL,
      Named("mean.coexp") = meanCor                                                                                                                                                                                                             
    );  
  }
  return List::create(
    Named("cor.coexp") = corCor,
    Named("mean.coexp") = meanCor
  );  
}

//' Calculate the Correlation of Coexpression and Mean Sign-Aware Coexpression
//'
//' Both of the coexpression statistics are calculated using all pairwise 
//' coexpression values in both the \emph{discovery} and \emph{test} datasets
//' respectively. For the other statistics, it makes sense to calculate the 
//' properties for the discovery network in advance to reduce calculation time
//' and memory. For these statistics this strategy doesn't make sense.
//'
//' @param pCoexpD,pCoexpT SEXP containers for the pointers to the coexpression 
//'  matrices for the \emph{discovery} and \emph{test} networks respectively.
//' @param discIndices,testIndices indices of the network subset in
//'   the \emph{discovery} and \emph{test} networks respectively.
//'
//' @return
//'   A vector containing:
//'   \enumerate{
//'     \item{\emph{cor.coexp}:}{
//'       The correlation between the subset coexpression for both networks.
//'     }
//'     \item{\emph{mean.coexp}:}{
//'       The mean correlation density of the network subset.
//'     }
//'   }
//'   
//' @references
//'   \enumerate{
//'     \item{
//'       Langfelder, P., Luo, R., Oldham, M. C. & Horvath, S. \emph{Is my 
//'       network module preserved and reproducible?} PLoS Comput. Biol. 
//'       \strong{7}, e1001057 (2011). 
//'     }
//'  }
//' @rdname CoexpStats-cpp
// [[Rcpp::export]]
List CoexpStats(
  SEXP pCoexpD, IntegerVector discIndices,
  SEXP pCoexpT, IntegerVector testIndices
) {
  XPtr<BigMatrix> xpCoexpD(pCoexpD);
  XPtr<BigMatrix> xpCoexpT(pCoexpT);

  if (xpCoexpD->ncol() != xpCoexpD->nrow()) {
    throw Rcpp::exception("provided 'pCoexpD' matrix is not square!");
  }
  if (xpCoexpT->ncol() != xpCoexpT->nrow()) {
    throw Rcpp::exception("provided 'pCoexpT' matrix is not square!");
  }

  if (discIndices.size() != testIndices.size()) {
    throw Rcpp::exception("'discIndices' and 'testIndices' must be the same length!");
  }

  // Make sure we're not indexing out of range.
  if (is_true(any(discIndices <= 0)) ||
      is_true(any(discIndices > xpCoexpD->ncol()))) {
    throw std::out_of_range(
      "Some of the requested indices for the discovery network are out of range!"
    );
  }
  if (is_true(any(testIndices <= 0)) ||
     is_true(any(testIndices > xpCoexpT->ncol()))) {
    throw std::out_of_range(
      "Some of the requested indices for the test network are out of range!"
    );
  }

  //  Dispatch function for all types of big.matrix.
  unsigned short typeD = xpCoexpD->matrix_type();
  unsigned short typeT = xpCoexpT->matrix_type();
  
  if (typeD == 1) {
    if (typeT == 1) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<char>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<char>(*xpCoexpT), testIndices
      );
    } else if (typeT == 2) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<char>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<short>(*xpCoexpT), testIndices
      );
    } else if (typeT == 4) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<char>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<int>(*xpCoexpT), testIndices
      );
    } else if (typeT == 8) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<char>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<double>(*xpCoexpT), testIndices
      );
    } else {
      /* We should never get here, unless the underlying implementation of
      bigmemory changes */
      throw Rcpp::exception("Undefined type for pCoexpT 'big.matrix' object");
    }
  } else if (typeD == 2) {
    if (typeT == 1) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<short>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<char>(*xpCoexpT), testIndices
      );
    } else if (typeT == 2) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<short>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<short>(*xpCoexpT), testIndices
      );
    } else if (typeT == 4) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<short>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<int>(*xpCoexpT), testIndices
      );
    } else if (typeT == 8) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<short>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<double>(*xpCoexpT), testIndices
      );
    } else {
      /* We should never get here, unless the underlying implementation of
      bigmemory changes */
      throw Rcpp::exception("Undefined type for pCoexpT 'big.matrix' object");
    }
  } else if (typeD == 4) {
    if (typeT == 1) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<int>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<char>(*xpCoexpT), testIndices
      );
    } else if (typeT == 2) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<int>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<short>(*xpCoexpT), testIndices
      );
    } else if (typeT == 4) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<int>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<int>(*xpCoexpT), testIndices
      );
    } else if (typeT == 8) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<int>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<double>(*xpCoexpT), testIndices
      );
    } else {
      /* We should never get here, unless the underlying implementation of
      bigmemory changes */
      throw Rcpp::exception("Undefined type for pCoexpT 'big.matrix' object");
    }
  } else if (typeD == 8) {
    if (typeT == 1) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<double>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<char>(*xpCoexpT), testIndices
      );
    } else if (typeT == 2) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<double>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<short>(*xpCoexpT), testIndices
      );
    } else if (typeT == 4) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<double>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<int>(*xpCoexpT), testIndices
      );
    } else if (typeT == 8) {
      return CoexpStats(
        xpCoexpD, MatrixAccessor<double>(*xpCoexpD), discIndices,
        xpCoexpT, MatrixAccessor<double>(*xpCoexpT), testIndices
      );
    } else {
      /* We should never get here, unless the underlying implementation of
      bigmemory changes */
      throw Rcpp::exception("Undefined type for pCoexpT 'big.matrix' object");
    }
  } else {
    /* We should never get here, unless the underlying implementation of
    bigmemory changes */
    throw Rcpp::exception("Undefined type for pCoexpD 'big.matrix' object");
  }
}

