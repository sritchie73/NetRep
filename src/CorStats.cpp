#define ARMA_USE_LAPACK
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
#include <bigmemory/MatrixAccessor.hpp>

// Implementation of CorStats
//
// @param xpCorD,xpCorT external pointers to the 'bigMatrix' objects for 
//  the discovery and test correlation structure matrices respectively.
// @param matCorD,matCorT Matrix Accessors for the 'bigMatrix' objects for
//  the discovery and test correlation structure matrices respectively.
// @param dIdx,tIdx indices for the variables composign the module of interest 
//  in the discovery and test datasets respectively.
//
// @return A list containing the
//    - a flattened vector of the module's correlation structure in the
//     discovery dataset (cor.discovery).
//    - a flattened vector of the module's correlation structure in the 
//      test dataset (cor.test).
//    - sign aware mean of the correlation (mean.cor)
template <typename S, typename T>
List CorStats(
  XPtr<BigMatrix> xpCorD, MatrixAccessor<S> matCorD, IntegerVector dIdx,
  XPtr<BigMatrix> xpCorT, MatrixAccessor<T> matCorT, IntegerVector tIdx
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
      corD.at(vi) = (double)matCorD[dIdx[jj] - 1][dIdx[ii] - 1]; 
      corT.at(vi) = (double)matCorT[tIdx[jj] - 1][tIdx[ii] - 1]; 
      vi++;
    }   
  }

  mat meanCor = mean(sign(corD) % corT, 1); 

  return List::create(
    Named("cor.discovery") = corD,
    Named("cor.test") = corT,
    Named("mean.cor") = meanCor
  );  
}

//' Calculate the correlation based statistics
//'
//' Both of the correlation statistics are calculated using all pairwise 
//' correlation values in both the \emph{discovery} and \emph{test} datasets.
//' For the other statistics, it makes sense to calculate the 
//' properties for the discovery network in advance to reduce calculation time
//' and memory. However, for the correaltion statistics this strategy doesn't 
//' make sense since we'd have to store a huge component of the discovery 
//' correlation structure.
//'
//' @param pCorD,pCorT SEXP containers for the pointers to the correlation 
//'  structure matrices for the \emph{discovery} and \emph{test} networks.
//' @param discIndices,testIndices indices of the network subset in
//'   the \emph{discovery} and \emph{test} networks respectively.
//'
//' @return
//'   A vector containing:
//'   \enumerate{
//'     \item{\emph{cor.discovery}:}{
//'       A flattened vector of the module's correlation structure in the 
//'       \emph{discovery} dataset.
//'     }
//'     \item{\emph{cor.test}:}{
//'       A flattened vector of the module's correlation structure in the 
//'       \emph{test} dataset.
//'     }
//'     \item{\emph{mean.cor}:}{
//'       The mean sign-aware correlation density of the network module.
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
//' @rdname CorStats-cpp
// [[Rcpp::export]]
List CorStats(
  SEXP pCorD, IntegerVector discIndices,
  SEXP pCorT, IntegerVector testIndices
) {
  XPtr<BigMatrix> xpCorD(pCorD);
  XPtr<BigMatrix> xpCorT(pCorT);

  if (xpCorD->ncol() != xpCorD->nrow()) {
    throw Rcpp::exception("provided 'pCorD' matrix is not square!");
  }
  if (xpCorT->ncol() != xpCorT->nrow()) {
    throw Rcpp::exception("provided 'pCorT' matrix is not square!");
  }

  if (discIndices.size() != testIndices.size()) {
    throw Rcpp::exception("'discIndices' and 'testIndices' must be the same length!");
  }

  // Make sure we're not indexing out of range.
  if (is_true(any(discIndices <= 0)) ||
      is_true(any(discIndices > xpCorD->ncol()))) {
    throw std::out_of_range(
      "Some of the requested indices for the discovery network are out of range!"
    );
  }
  if (is_true(any(testIndices <= 0)) ||
     is_true(any(testIndices > xpCorT->ncol()))) {
    throw std::out_of_range(
      "Some of the requested indices for the test network are out of range!"
    );
  }

  //  Dispatch function for all types of big.matrix.
  unsigned short typeD = xpCorD->matrix_type();
  unsigned short typeT = xpCorT->matrix_type();
  
  if (typeD == 1) {
    if (typeT == 1) {
      return CorStats(
        xpCorD, MatrixAccessor<char>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<char>(*xpCorT), testIndices
      );
    } else if (typeT == 2) {
      return CorStats(
        xpCorD, MatrixAccessor<char>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<short>(*xpCorT), testIndices
      );
    } else if (typeT == 4) {
      return CorStats(
        xpCorD, MatrixAccessor<char>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<int>(*xpCorT), testIndices
      );
    } else if (typeT == 6) {
      return CorStats(
        xpCorD, MatrixAccessor<char>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<float>(*xpCorT), testIndices
      );
    }else if (typeT == 8) {
      return CorStats(
        xpCorD, MatrixAccessor<char>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<double>(*xpCorT), testIndices
      );
    } else {
      /* We should never get here, unless the underlying implementation of
      bigmemory changes */
      throw Rcpp::exception("Undefined type for pCorT 'bigMatrix' object");
    }
  } else if (typeD == 2) {
    if (typeT == 1) {
      return CorStats(
        xpCorD, MatrixAccessor<short>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<char>(*xpCorT), testIndices
      );
    } else if (typeT == 2) {
      return CorStats(
        xpCorD, MatrixAccessor<short>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<short>(*xpCorT), testIndices
      );
    } else if (typeT == 4) {
      return CorStats(
        xpCorD, MatrixAccessor<short>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<int>(*xpCorT), testIndices
      );
    } else if (typeT == 6) {
      return CorStats(
        xpCorD, MatrixAccessor<short>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<float>(*xpCorT), testIndices
      );
    } else if (typeT == 8) {
      return CorStats(
        xpCorD, MatrixAccessor<short>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<double>(*xpCorT), testIndices
      );
    } else {
      /* We should never get here, unless the underlying implementation of
      bigmemory changes */
      throw Rcpp::exception("Undefined type for 'pCorT' 'bigMatrix' object");
    }
  } else if (typeD == 4) {
    if (typeT == 1) {
      return CorStats(
        xpCorD, MatrixAccessor<int>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<char>(*xpCorT), testIndices
      );
    } else if (typeT == 2) {
      return CorStats(
        xpCorD, MatrixAccessor<int>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<short>(*xpCorT), testIndices
      );
    } else if (typeT == 4) {
      return CorStats(
        xpCorD, MatrixAccessor<int>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<int>(*xpCorT), testIndices
      );
    } else if (typeT == 6) {
      return CorStats(
        xpCorD, MatrixAccessor<int>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<float>(*xpCorT), testIndices
      );
    } else if (typeT == 8) {
      return CorStats(
        xpCorD, MatrixAccessor<int>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<double>(*xpCorT), testIndices
      );
    } else {
      /* We should never get here, unless the underlying implementation of
      bigmemory changes */
      throw Rcpp::exception("Undefined type for 'pCorT' 'bigMatrix' object");
    }
  } else if (typeD == 6) {
    if (typeT == 1) {
      return CorStats(
        xpCorD, MatrixAccessor<float>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<char>(*xpCorT), testIndices
      );
    } else if (typeT == 2) {
      return CorStats(
        xpCorD, MatrixAccessor<float>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<short>(*xpCorT), testIndices
      );
    } else if (typeT == 4) {
      return CorStats(
        xpCorD, MatrixAccessor<float>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<int>(*xpCorT), testIndices
      );
    } else if (typeT == 6) {
      return CorStats(
        xpCorD, MatrixAccessor<float>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<float>(*xpCorT), testIndices
      );
    } else if (typeT == 8) {
      return CorStats(
        xpCorD, MatrixAccessor<float>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<double>(*xpCorT), testIndices
      );
    } else {
      /* We should never get here, unless the underlying implementation of
      bigmemory changes */
      throw Rcpp::exception("Undefined type for 'pCorT' 'bigMatrix' object");
    }
  } else if (typeD == 8) {
    if (typeT == 1) {
      return CorStats(
        xpCorD, MatrixAccessor<double>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<char>(*xpCorT), testIndices
      );
    } else if (typeT == 2) {
      return CorStats(
        xpCorD, MatrixAccessor<double>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<short>(*xpCorT), testIndices
      );
    } else if (typeT == 4) {
      return CorStats(
        xpCorD, MatrixAccessor<double>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<int>(*xpCorT), testIndices
      );
    } else if (typeT == 6) {
      return CorStats(
        xpCorD, MatrixAccessor<double>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<float>(*xpCorT), testIndices
      );
    } else if (typeT == 8) {
      return CorStats(
        xpCorD, MatrixAccessor<double>(*xpCorD), discIndices,
        xpCorT, MatrixAccessor<double>(*xpCorT), testIndices
      );
    } else {
      /* We should never get here, unless the underlying implementation of
      bigmemory changes */
      throw Rcpp::exception("Undefined type for 'pCorT' 'bigMatrix' object");
    }
  } else {
    /* We should never get here, unless the underlying implementation of
    bigmemory changes */
    throw Rcpp::exception("Undefined type for 'pCorD' 'bigMatrix' object");
  }
}

