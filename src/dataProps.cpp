#define ARMA_USE_LAPACK
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(BH, bigmemory, RcppArmadillo)]]
#include <bigmemory/BigMatrix.h>

/* Implementation of DataProps
 *
 * @param dat the armadillo compatible gene expression matrix
 * @param subsetIndices indices of the network subset of interest.
 * @return
 *    A List containing:
 *     - A vector summarising the module (moduleSummary).
 *     - The Module Membership of each variable (MM).
 *     - The proportion of module variance explained by the module summary vector 
 *       (pve).
 */
template <typename T>
List DataProps(const Mat<T>& dat, IntegerVector subsetIndices) {
  Mat<T> U, V;
  Col<T> S;
  uvec subsetCols = sort(as<uvec>(subsetIndices) - 1);
  
  // Get the summary profile for the network subset from the SVD.
  bool success = svd_econ(U, S, V, dat.cols(subsetCols), "left", "dc");
  if (!success) {
    Function warning("warning");
    warning("SVD failed to converge, does your data contain missing or"
            " infinite values?");
    return List::create(
        Named("SP") = NA_REAL,
        Named("MM") = NA_REAL,
        Named("pve") = NA_REAL
      );
  }
  Mat<T> summary(U.col(0));

  // Flip the sign of the summary profile so that the eigenvector is 
  // positively correlated with the average scaled value of the underlying
  // data for the network module.
  Mat<T> ap = cor(mean(dat.cols(subsetCols), 1), summary);
  if (ap(0,0) < 0) {
    summary *= -1; 
  }
  
  // We want the correlation between each variable (node) in the underlying
  // data and the summary profile for that network subset.
  Mat<T> p = cor(summary, dat.cols(subsetCols));
  Mat<T> MM(p);
  
  // To make sure the resulting KIM vector is in the correct order,
  // order the results to match the original ordering of subsetIndices.
  Function rank("rank"); // Rank only works on R objects like IntegerVector.
  uvec idxRank = as<uvec>(rank(subsetIndices)) - 1;

  Col<T> oMM = MM(idxRank);

  // The proportion of variance explained is the sum of the squared 
  // correlation between the network module summary profile, and each of the 
  // variables in the data that correspond to nodes in the network subset.
  Mat<T> pve(mean(square(p), 1));
  
  return List::create(
    Named("moduleSummary") = summary,
    Named("MM") = oMM,
    Named("pve") = pve
  );
}

//' Get the module summary vector and the proportion of variance it explains
//' 
//' @param pDat SEXP container for the pointer to a scaled version of the 
//'   data matrix used to construct the network.
//' @param subsetIndices indices of the network subset of interest in 
//'   \code{pDat}.
//' 
//' @return
//'  A list containing:
//'  \enumerate{
//'   \item{\emph{"moduleSummary"}:}{
//'     The module summary profile (see details).
//'   }
//'   \item{\emph{"MM"}:}{
//'     The Module Membership of each node (see details).
//'   }
//'   \item{\emph{"pve"}:}{
//'     The proportion of the variance explained by the module's summary 
//'      profile (see details).
//'   }
//'  }
//'  
//' @references
//'  \enumerate{
//'     \item{
//'       Langfelder, P., Luo, R., Oldham, M. C. & Horvath, S. \emph{Is my
//'       network module preserved and reproducible?} PLoS Comput. Biol. 
//'       \strong{7}, e1001057 (2011). 
//'     }
//'  }
//'  
//' @details
//'  First, the module summary profile (SP) is calculated as the first 
//'  eigenvector of a principal component analysis of the variables composing 
//'  the module of interest. The orientation of the eigenvector is modified so 
//'  that its sign is in the same direction as the module on average. I.e. for 
//'  gene expression data this is the "module eigengene" \emph{(1)}.
//'  
//'  The Module Membership (MM) is the correlation between each variable 
//'  composing the module and the module's summary profile.
//'  
//'  The proportion of module variance explained by the summary profile (pve) 
//'  is quantified as the average square of the Module Membership \emph{(1)}.
//' 
//' @import RcppArmadillo
//' @rdname dataProps-cpp
//'  
// [[Rcpp::export]]
List DataProps(SEXP pDat, IntegerVector subsetIndices) {
  XPtr<BigMatrix> xpDat(pDat);
  
  // Make sure we're not indexing out of range.
  if (is_true(any(subsetIndices <= 0)) || 
      is_true(any(subsetIndices > xpDat->ncol()))) {
    throw std::out_of_range("Some of requested indices are outside of range!");
  }
  
  // Dispatch function for all types of big.matrix.
  unsigned short type = xpDat->matrix_type();
  if (type == 6) {
    return DataProps(
      arma::Mat<float>((float *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false),
      subsetIndices
    );
  } else if (type == 8) {
    return DataProps(
      arma::Mat<double>((double *)xpDat->matrix(), xpDat->nrow(), xpDat->ncol(), false),
      subsetIndices
    );
  } else {
    throw Rcpp::exception(
      "SVD can only be calculated on a 'bigMatrix' whose underlying type is"
      "'double' or 'float'."
    );
  }
}