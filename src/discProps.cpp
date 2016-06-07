#include "utils.h"
#include "netStats.h"

//' Calculate the intermediate network properties in the discovery dataset
//' 
//' These properties are need at every permutation: so they will be computed 
//' once.
//' 
//' @details
//' \subsection{Input expectations:}{
//'   Note that this function expects all inputs to be sensible, as checked by
//'   the R function 'checkUserInput' and processed by 'modulePreservation'. 
//'   
//'   These requirements are:
//'   \itemize{
//'   \item{The ordering of node names across 'dData', 'dCorr', and 'dNet' is
//'         consistent.}
//'   \item{The columns of 'dData' are the nodes.}
//'   \item{'dCorr' and 'dNet'  are square matrices, and their rownames are 
//'         identical to their column names.}
//'   \item{'moduleAssigments' is a named character vector, where the names
//'         represent node labels found in the discovery dataset (e.g. 'dNet').}
//'   }
//' }
//' 
//' @param dData data matrix from the \emph{discovery} dataset.
//' @param dCorr matrix of correlation coefficients between all pairs of 
//'   variables/nodes in the \emph{discovery} dataset.
//' @param dNet adjacency matrix of network edge weights between all pairs of 
//'   nodes in the \emph{discovery} dataset.
//' @param tNodeNames a character vector of node names in the test dataset
//' @param moduleAssignments a named character vector containing the module 
//'   each node belongs to in the discovery dataset. 
//' @param modules a character vector of modules for which to calculate the 
//'   module preservation statistics.
//' 
//' @return a list containing three lists: a list of weighted degree vectors,
//'   a list of correlation coefficient vectors, and a list of node 
//'   contribution vectors. There is one vector for each module in each list.
//'   
// [[Rcpp::export]]
Rcpp::List IntermediateProperties (
    Rcpp::NumericMatrix dData, Rcpp::NumericMatrix dCorr, Rcpp::NumericMatrix dNet,
    Rcpp::CharacterVector tNodeNames, Rcpp::CharacterVector moduleAssignments, 
    Rcpp::CharacterVector modules
) {
  // First, we need to create pointers to the memory holding each
  // NumericMatrix that can be recognised by the armadillo library.
  arma::mat dDataPtr = arma::mat(dData.begin(), dData.nrow(), dData.ncol(), false, true);
  arma::mat dCorrPtr = arma::mat(dCorr.begin(), dCorr.nrow(), dCorr.ncol(), false, true);
  arma::mat dNetPtr = arma::mat(dNet.begin(), dNet.nrow(), dNet.ncol(), false, true);
  
  // Next we will scale the matrix data
  arma::mat dDataScaled = Scale(dDataPtr);

  R_CheckUserInterrupt(); 
  
  // convert the colnames / rownames to C++ equivalents
  const std::vector<std::string> dNames (Rcpp::as<std::vector<std::string>>(colnames(dNet)));
  const std::vector<std::string> tNames (Rcpp::as<std::vector<std::string>>(tNodeNames));
  
  /* Next, we need to create three mappings:
  *  - From node IDs to indices in the discovery dataset.
  *  - From modules to all node IDs.
  *  - From modules to just node IDs present in the test dataset.
  */
  const namemap dIdxMap = MakeIdxMap(dNames);
  const stringmap modNodeMap = MakeModMap(moduleAssignments);
  const namemap tIdxMap = MakeIdxMap(tNames);
  const stringmap modNodePresentMap = MakeModMap(moduleAssignments, tIdxMap);
  
  // What modules do we actually want to analyse?
  const std::vector<std::string> mods (Rcpp::as<std::vector<std::string>>(modules));
  
  // We only need to iterate through modules which have nodes in the test 
  // dataset
  std::vector<std::string> modsPresent;
  for (auto it = mods.begin(); it != mods.end(); ++it) {
    if (modNodePresentMap.count(*it) > 0) {
      modsPresent.push_back(*it);
    }
  }
  
  R_CheckUserInterrupt(); 
  
  Rcpp::List degree;
  Rcpp::List corr;
  Rcpp::List contribution;
  
  // Calculate the network properties in the discovery dataset.
  std::string mod;
  arma::uvec dIdx, dRank;
  arma::vec dSP, dWD, dCV, dNC; 
  for (auto mi = modsPresent.begin(); mi != modsPresent.end(); ++mi) {
    // Get the node indices in the discovery dataset for this module
    mod = *mi;
    dIdx = GetNodeIdx(mod, modNodePresentMap, dIdxMap);
    R_CheckUserInterrupt(); 
    
    // Calculate the network properties and insert into their storage containers
    dCV = CorrVector(dCorrPtr, dIdx);
    
    R_CheckUserInterrupt(); 
    
    dRank = sortNodes(dIdx); // Sort nodes for sequential memory access
    dWD = WeightedDegree(dNetPtr, dIdx)(dRank);
    R_CheckUserInterrupt(); 
    dSP = SummaryProfile(dDataScaled, dIdx);
    R_CheckUserInterrupt(); 
    dNC = NodeContribution(dDataScaled, dIdx, dSP)(dRank);
    R_CheckUserInterrupt(); 
    
    // Cast to R-vectors and add to results lists
    corr.push_back(Rcpp::NumericVector(dCV.begin(), dCV.end()));
    degree.push_back(Rcpp::NumericVector(dWD.begin(), dWD.end()));
    contribution.push_back(Rcpp::NumericVector(dNC.begin(), dNC.end()));
  }
  degree.names() = modsPresent;
  corr.names() = modsPresent;
  contribution.names() = modsPresent;
  
  return Rcpp::List::create(
    Rcpp::Named("degree") = degree,
    Rcpp::Named("corr") = corr,
    Rcpp::Named("contribution") = contribution
  );
}

//' Calculate the intermediate network properties in the discovery dataset
//' 
//' These properties are need at every permutation: so they will be computed 
//' once.
//' 
//' @details
//' \subsection{Input expectations:}{
//'   Note that this function expects all inputs to be sensible, as checked by
//'   the R function 'checkUserInput' and processed by 'modulePreservation'. 
//'   
//'   These requirements are:
//'   \itemize{
//'   \item{The ordering of node names across 'dCorr' and 'dNet' is
//'         consistent.}
//'   \item{'dCorr' and 'dNet'  are square matrices, and their rownames are 
//'         identical to their column names.}
//'   \item{'moduleAssigments' is a named character vector, where the names
//'         represent node labels found in the discovery dataset (e.g. 'dNet').}
//'   }
//' }
//' 
//' @param dCorr matrix of correlation coefficients between all pairs of 
//'   variables/nodes in the \emph{discovery} dataset.
//' @param dNet adjacency matrix of network edge weights between all pairs of 
//'   nodes in the \emph{discovery} dataset.
//' @param tNodeNames a character vector of node names in the test dataset
//' @param moduleAssignments a named character vector containing the module 
//'   each node belongs to in the discovery dataset. 
//' @param modules a character vector of modules for which to calculate the 
//'   module preservation statistics.
//' 
//' @return a list containing two lists: a list of weighted degree vectors,
//'   and a list of correlation coefficient vectors. There is one vector for 
//'   each module in each list.
//'   
// [[Rcpp::export]]
Rcpp::List IntermediatePropertiesNoData (
    Rcpp::NumericMatrix dCorr, Rcpp::NumericMatrix dNet,
    Rcpp::CharacterVector tNodeNames, Rcpp::CharacterVector moduleAssignments, 
    Rcpp::CharacterVector modules
) {
  // First, we need to create pointers to the memory holding each
  // NumericMatrix that can be recognised by the armadillo library.
  arma::mat dCorrPtr = arma::mat(dCorr.begin(), dCorr.nrow(), dCorr.ncol(), false, true);
  arma::mat dNetPtr = arma::mat(dNet.begin(), dNet.nrow(), dNet.ncol(), false, true);
  
  R_CheckUserInterrupt(); 
  
  // convert the colnames / rownames to C++ equivalents
  const std::vector<std::string> dNames (Rcpp::as<std::vector<std::string>>(colnames(dNet)));
  const std::vector<std::string> tNames (Rcpp::as<std::vector<std::string>>(tNodeNames));
  
  /* Next, we need to create three mappings:
  *  - From node IDs to indices in the discovery dataset.
  *  - From modules to all node IDs.
  *  - From modules to just node IDs present in the test dataset.
  */
  const namemap dIdxMap = MakeIdxMap(dNames);
  const stringmap modNodeMap = MakeModMap(moduleAssignments);
  const namemap tIdxMap = MakeIdxMap(tNames);
  const stringmap modNodePresentMap = MakeModMap(moduleAssignments, tIdxMap);
  
  // What modules do we actually want to analyse?
  const std::vector<std::string> mods (Rcpp::as<std::vector<std::string>>(modules));

  // We only need to iterate through modules which have nodes in the test 
  // dataset
  std::vector<std::string> modsPresent;
  for (auto it = mods.begin(); it != mods.end(); ++it) {
    if (modNodePresentMap.count(*it) > 0) {
      modsPresent.push_back(*it);
    }
  }
  
  R_CheckUserInterrupt(); 
  
  Rcpp::List degree;
  Rcpp::List corr;

  // Calculate the network properties in the discovery dataset.
  std::string mod;
  arma::uvec dIdx, dRank;
  arma::vec dWD, dCV; 
  for (auto mi = modsPresent.begin(); mi != modsPresent.end(); ++mi) {
    // Get the node indices in the discovery dataset for this module
    mod = *mi;
    dIdx = GetNodeIdx(mod, modNodePresentMap, dIdxMap);
    R_CheckUserInterrupt(); 
    
    // Calculate the network properties and insert into their storage containers
    dCV = CorrVector(dCorrPtr, dIdx);
    
    R_CheckUserInterrupt(); 
    
    dRank = sortNodes(dIdx); // Sort nodes for sequential memory access
    dWD = WeightedDegree(dNetPtr, dIdx)(dRank);
    R_CheckUserInterrupt(); 
    
    // Cast to R-vectors and add to results lists
    corr.push_back(Rcpp::NumericVector(dCV.begin(), dCV.end()));
    degree.push_back(Rcpp::NumericVector(dWD.begin(), dWD.end()));
  }
  degree.names() = modsPresent;
  corr.names() = modsPresent;

  return Rcpp::List::create(
    Rcpp::Named("degree") = degree,
    Rcpp::Named("corr") = corr
  );
}
