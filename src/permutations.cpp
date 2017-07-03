#include "utils.h"
#include "netStats.h"
#include "thread-utils.h"

/* Generate null-distribution observations for the module preservation statistics
 * 
 * Fills out the corresponding slices of the provided 'nulls' cube based on the
 * number of permutations requested, and the start index.
 * 
 * @param tDataAddr memory address of the (scaled) test data matrix.
 * @param tCorrAddr memory address of the test correlation matrix.
 * @param tNetAddr memory address of the test network matrix.
 * @param nSamples number of samples in the test dataset.
 * @param nNodes number of nodes in the test network.
 * @param addrWD unordered_map containing the memory addresses of the weighted 
 *   degree vectors for all modules in the discovery dataset.
 * @param addrNC unordered_map containing the memory addresses of the node 
 *   contribution vectors for all modules in the discovery dataset.
 * @param addrCV unordered_map containing the memory addresses of the 
 *   correlation coefficient vectors for all modules in the discovery dataset.
 * @param mods vector of modules for which the module preservation statistics
 *   are being calculated for.
 * @param modNodeMap mapping of module labels to node IDs.
 * @param modIdxMap mapping of module labels to 'nulls' cube row indices.
 * @param nullIdx a vector of node IDs to be shuffled in the permutation 
 *  procedure.
 * @param nullMap mapping of node IDs to indices in 'nullIdx'.
 * @param nullsAddr memory address of the cube to store the results in
 * @param totalPerm total number of permutations.
 * @param nPerm number of permutations for this thread to calculate.
 * @param start slice index to start at when filling in the 'nulls' cube.
 * @param progressAddr memory address of the vector to fill in number of 
 *  permutations completed for this thread.
 * @param nThreads total number of threads executing.
 * @param thread the number of the thread.
 * @param interrupted variable on the heap checking whether the user has asked
 *   to cancel the computation
 */
void calculateNulls(
  double * tDataAddr, double * tCorrAddr, double * tNetAddr, 
  unsigned int nSamples, unsigned int nNodes, addrmap& addrWD, addrmap& addrNC, 
  addrmap& addrCV, const std::vector<std::string> mods, 
  const stringmap modNodeMap, const namemap modIdxMap, arma::uvec nullIdx,
  namemap nullMap, double * nullsAddr, unsigned int totalPerm, 
  unsigned int nPerm, unsigned int start, unsigned int * progressAddr,
  unsigned int nThreads, unsigned int thread, bool& interrupted
) {    
  /**
   * Note: the R API is single threaded, we *must not* access it
   * at all in this function or any functions it calls (i.e. netStats.cpp).
   **/
  
  // Tell this thread where the results cube and progress bar are located in 
  // memory:
  arma::cube nulls = arma::cube(nullsAddr, mods.size(), 7, totalPerm, false, true);
  arma::uvec progress = arma::uvec(progressAddr, nThreads, false, true);
  
  std::string mod;
  unsigned int modIdx, mNodes;
  arma::uvec tIdx, tRank;
  arma::vec tWD, tSP, tNC, tCV;
  for (unsigned int pp = start; pp < start + nPerm; ++pp) {
    nullIdx = arma::shuffle(nullIdx);
    for (auto mi = mods.begin(); mi != mods.end(); ++mi) {
      if (interrupted) return; 
      // What module are we analysing, and what index does it have internally?
      mod = *mi; 
      modIdx = modIdxMap.at(mod);
     
      // Get the node indices in the test dataset for this module
      tIdx = GetRandomIdx(mod, modNodeMap, nullIdx.memptr(), nullIdx.n_elem, nullMap);
      mNodes = tIdx.n_elem;
      
      // Now calculate required properties in the test dataset
      tCV = CorrVector(tCorrAddr, nNodes, tIdx.memptr(), mNodes);
      if (interrupted) return; 
      
      // Sort nodes indices for sequential memory access
      tRank = SortNodes(tIdx.memptr(), mNodes); 
      
      tWD = WeightedDegree(tNetAddr, nNodes, tIdx.memptr(), mNodes);
      tWD = tWD(tRank); //reorder results
      if (interrupted) return; 
      
      tSP = SummaryProfile(tDataAddr, nSamples, nNodes, tIdx.memptr(), mNodes);
      if (interrupted) return; 
      
      tNC = NodeContribution(tDataAddr, nSamples, nNodes, tIdx.memptr(),
                             mNodes, tSP.memptr());
      tNC = tNC(tRank); // reorder results
      if (interrupted) return; 
      
      // Calculate and store test statistics in the appropriate location in the 
      // results matrix
      nulls.at(modIdx, 0, pp) = AverageEdgeWeight(tWD.memptr(), tWD.n_elem);
      nulls.at(modIdx, 1, pp) = ModuleCoherence(tNC.memptr(), tNC.n_elem);
      nulls.at(modIdx, 2, pp) = Correlation(addrCV[mod], tCV.memptr(), tCV.n_elem);
      nulls.at(modIdx, 3, pp) = Correlation(addrWD[mod], tWD.memptr(), tWD.n_elem);
      nulls.at(modIdx, 4, pp) = Correlation(addrNC[mod], tNC.memptr(), tNC.n_elem);
      nulls.at(modIdx, 5, pp) = SignAwareMean(addrCV[mod], tCV.memptr(), tCV.n_elem);
      nulls.at(modIdx, 6, pp) = SignAwareMean(addrNC[mod], tNC.memptr(), tNC.n_elem);
    }
    progress[thread]++; 
  }
}

///' Multithreaded permutation procedure for module preservation statistics
///' 
///' @details
///' \subsection{Input expectations:}{
///'   Note that this function expects all inputs to be sensible, as checked by
///'   the R function 'checkUserInput' and processed by 'modulePreservation'. 
///'   
///'   These requirements are:
///'   \itemize{
///'   \item{The ordering of node names across 'tData', 'tCorr', and 'tNet' is
///'         consistent.}
///'   \item{The columns of 'tData' are the nodes.}
///'   \item{'tData' has been scaled by 'Scale'.}
///'   \item{'tCorr' and 'tNet' are square matrices, and their rownames are 
///'         identical to their column names.}
///'   \item{'moduleAssigments' is a named character vector, where the names
///'         represent node labels found in the discovery dataset (e.g. 'dNet').}
///'   \item{'nPermutations' is a single number, greater than 0.}
///'   \item{'nCores' is a single number, greater than 0. Note, this number must
///'         not be larger than the number of cores on your machine, or the 
///'         number of cores allocated to your job!}
///'   \item{'nullHypothesis' must be a character vector of length 1, containing
///'         either "overlap" or "all".}
///'   \item{'verbose' must be a logical vector of length 1 containing either 
///'         'TRUE' or 'FALSE'.}
///'   \item{'vCat' must be the function NetRep:::vCat.}
///'   }
///' }
///' 
///' @param discProps a list of intermediate properties calculated in the 
///'   discovery dataset by \code{\link{IntermediateProperties}}.
///' @param tData scaled data matrix from the \emph{test} dataset.
///' @param tCorr matrix of correlation coefficients between all pairs of 
///'   variables/nodes in the \emph{test} dataset.
///' @param tNet adjacency matrix of network edge weights between all pairs of 
///'   nodes in the \emph{test} dataset.
///' @param moduleAssignments a named character vector containing the module 
///'   each node belongs to in the discovery dataset. 
///' @param modules a character vector of modules for which to calculate the 
///'   module preservation statistics.
///' @param nPermutations the number of permutations from which to generate the
///'   null distributions for each statistic.
///' @param nCores the number of cores that the permutation procedure may use.
///' @param nullHypothesis either "overlap" or "all".
///' @param verbose if 'true', then progress messages are printed.
///' @param vCat the vCat function must be passed in so that it can be called 
///'  for output logging. 
///' 
///' @return a list containing a matrix of observed test statistics, and an
///'   array of null distribution observations.
///'   
///' @keywords internal
// [[Rcpp::export]]
Rcpp::List PermutationProcedure (
  Rcpp::List discProps, Rcpp::NumericMatrix tData, Rcpp::NumericMatrix tCorr, 
  Rcpp::NumericMatrix tNet, Rcpp::CharacterVector moduleAssignments, 
  Rcpp::CharacterVector modules, Rcpp::IntegerVector nPermutations, 
  Rcpp::IntegerVector nCores, Rcpp::CharacterVector nullHypothesis, 
  Rcpp::LogicalVector verbose, Rcpp::Function vCat
) {
  // convert the colnames / rownames to C++ equivalents
  const std::vector<std::string> dNames (Rcpp::as<std::vector<std::string>>(moduleAssignments.names()));
  const std::vector<std::string> tNames (Rcpp::as<std::vector<std::string>>(colnames(tNet)));
  unsigned int nSamples = tData.nrow();
  unsigned int nNodes = tData.ncol();
  
  // Define statistic names
  const std::vector<std::string> statnames = {
    "avg.weight", "coherence", "cor.cor", "cor.degree", "cor.contrib", 
    "avg.cor", "avg.contrib"
  };
  
  /* Next, we need to create three mappings:
   *  - From node IDs to indices in the test dataset.
   *  - From modules to all node IDs.
   *  - From modules to just node IDs present in the test dataset.
   */
  const namemap tIdxMap = MakeIdxMap(tNames);
  const stringmap modNodeMap = MakeModMap(moduleAssignments);
  const stringmap modNodePresentMap = MakeModMap(moduleAssignments, tIdxMap);
  
  // What modules do we actually want to analyse?
  const std::vector<std::string> mods (Rcpp::as<std::vector<std::string>>(modules));
  
  // Map module labels to row indices in our armadillo matrices/cubes
  const namemap modIdxMap = MakeIdxMap(mods);
  
  // We only need to iterate through modules which have nodes in the test 
  // dataset
  std::vector<std::string> modsPresent;
  for (auto it = mods.begin(); it != mods.end(); ++it) {
    if (modNodePresentMap.count(*it) > 0) {
      modsPresent.push_back(*it);
    }
  }

  // Typecast function options from R's vectors to appropriate C++ scalar 
  // equivalents
  std::string nullType = Rcpp::as<std::string>(nullHypothesis[0]);
  unsigned int nThreads = nCores[0];
  unsigned int nPerm = nPermutations[0];
  const bool verboseFlag = verbose[0];
  
  // Initialise results container for storing observed test statistics
  arma::mat obs (mods.size(), 7);
  obs.fill(NA_REAL);

  /* We need to convert each 'discProps' list to a mapping from each module
   * to the address in memory its corresponding property vector
   */

  // We need to do a lot of casting to get list elements in C++!
  Rcpp::List lWD = Rcpp::as<Rcpp::List>(discProps["degree"]);
  Rcpp::List lCV = Rcpp::as<Rcpp::List>(discProps["corr"]);
  Rcpp::List lNC = Rcpp::as<Rcpp::List>(discProps["contribution"]);
  Rcpp::NumericVector vWD, vCV, vNC; 
  
  addrmap addrWD, addrNC, addrCV; 
  std::string mod;
  for (auto mi = modsPresent.begin(); mi != modsPresent.end(); ++mi) {
    mod = *mi;
    
    // Extract the numeric vectors
    vWD = Rcpp::as<Rcpp::NumericVector>(lWD[mod]);
    vCV = Rcpp::as<Rcpp::NumericVector>(lCV[mod]);
    vNC = Rcpp::as<Rcpp::NumericVector>(lNC[mod]);
    
    addrWD[mod] = vWD.begin();
    addrCV[mod] = vCV.begin();
    addrNC[mod] = vNC.begin();
  }
  
  R_CheckUserInterrupt();
  
  // Now calculate the observed test statistics
  vCat(verbose, 1, "Calculating observed test statistics...");
  unsigned int modIdx, mNodes;
  arma::uvec tIdx, tRank;
  arma::vec tCV, tWD, tSP, tNC;
  for (auto mi = modsPresent.begin(); mi != modsPresent.end(); ++mi) {
    // What module are we analysing, and what index does it have internally?
    mod = *mi;
    modIdx = modIdxMap.at(mod);
    
    // Get the node indices in the test dataset for this module
    tIdx = GetNodeIdx(mod, modNodePresentMap, tIdxMap);
    mNodes = tIdx.n_elem;
    
    // Now calculate required properties in the test dataset
    tCV = CorrVector(tCorr.begin(), nNodes, tIdx.memptr(), mNodes);
    R_CheckUserInterrupt();
    
    // Sort node indices for sequential memory access
    tRank = SortNodes(tIdx.memptr(), mNodes); 
    
    tWD = WeightedDegree(tNet.begin(), nNodes, tIdx.memptr(), mNodes);
    tWD = tWD(tRank); // reorder results
    R_CheckUserInterrupt(); 
    
    tSP = SummaryProfile(tData.begin(), nSamples, nNodes, tIdx.memptr(), 
                         mNodes);
    R_CheckUserInterrupt(); 
    
    tNC = NodeContribution(tData.begin(), nSamples, nNodes, 
                           tIdx.memptr(), mNodes, tSP.memptr());
    tNC = tNC(tRank); // reorder results
    R_CheckUserInterrupt(); 
    
    /* Calculate and store test statistics in the appropriate location in the 
    * results matrix
    */
    obs(modIdx, 0) = AverageEdgeWeight(tWD.memptr(), tWD.n_elem);
    obs(modIdx, 1) = ModuleCoherence(tNC.memptr(), tNC.n_elem);
    obs(modIdx, 2) = Correlation(addrCV[mod], tCV.memptr(), tCV.n_elem);
    obs(modIdx, 3) = Correlation(addrWD[mod], tWD.memptr(), tWD.n_elem);
    obs(modIdx, 4) = Correlation(addrNC[mod], tNC.memptr(), tNC.n_elem);
    obs(modIdx, 5) = SignAwareMean(addrCV[mod], tCV.memptr(), tCV.n_elem);
    obs(modIdx, 6) = SignAwareMean(addrNC[mod], tNC.memptr(), tNC.n_elem);
  }
  
  // Just return observed statistics if no permutations requested
  if (nPerm == 0) {
    // Convert any NaNs or Infinites to NA_REALs
    obs.elem(arma::find_nonfinite(obs)).fill(NA_REAL);
    
    // Convert matrix of observed test statistics into an R object before
    // returning
    Rcpp::NumericMatrix observed (obs.n_rows, obs.n_cols, obs.begin());
    colnames(observed) = Rcpp::CharacterVector(statnames.begin(), statnames.end());
    rownames(observed) = modules;
    
    return Rcpp::List::create(Rcpp::Named("observed") = observed);
  }

  // If there are permutations requested, proceed.
  
  // Initialise results container for storing the null distributions
  arma::cube nulls (mods.size(), 7, nPerm); 
  nulls.fill(NA_REAL);
  
  /* For the permutation procedure, we need to shuffle a vector of *valid*
  * indices in the test network: if the null hypothesis is "overlap" (the
  * default) then only nodes that are present in both the discovery and test
  * datasets are used to generate the null distributions.
  * 
  * So we need:
  *  - A *vector* of indices in the test dataset that can be shuffled
  *  - A mapping from the valid node IDs to their indices in the vector to 
  *    be shuffled
  */
  arma::uvec nullIdx;
  namemap nullMap;
  if (nullType == "overlap") {
    nullMap = MakeNullMap(dNames, tIdxMap, nullIdx);
  } else { // otherwise take all nodes
    nullMap = MakeNullMap(tNames, tIdxMap, nullIdx);
  }
  R_CheckUserInterrupt(); 
  
  if (nThreads == 1) {
    vCat(verbose, 1, "Generating null distributions from", nPerm, 
         "permutations using", nThreads, "thread...");
  } else {
    vCat(verbose, 1, "Generating null distributions from", nPerm, 
         "permutations using", nThreads, "threads...");
  }

  // Create nThreads to run the permutation procedure in parallel
  std::thread *tt = new std::thread[nThreads];

  // Determine the number of permutations for each thread
  arma::uvec chunkPerms (nThreads, arma::fill::zeros);
  for (unsigned int ii = 0; ii < nThreads; ++ii) {
    chunkPerms.at(ii) = nPerm / nThreads;
  }
  // When the number of permtutations cannot be distributed evenly, spread the
  // remainder out across threads.
  for (unsigned int ii = 0; ii < nPerm % nThreads; ++ii) {
    chunkPerms.at(ii)++;
  }

  // Set up the slice indices each thread should start at in the results cube
  arma::uvec startIdx (nThreads, arma::fill::zeros);
  for (unsigned int ii = 0; ii < nThreads; ++ii) {
    for (unsigned int jj = 0; jj < ii; ++jj) {
      startIdx.at(ii) += chunkPerms.at(jj);
    }
  }

  // Set up the progress bar
  arma::uvec progress (nThreads, arma::fill::zeros);

  // This variable stored on the heap will get modified to be 'true' of ^C is
  // entered in the R terminal (see 'MonitorProgress'). Each thread will check
  // this, then abort if it is 'true'
  bool interrupted = false;

  // Spawn the threads
  for (unsigned int ii = 0; ii < nThreads; ++ii) {
    tt[ii] = std::thread(
      calculateNulls, tData.begin(), tCorr.begin(),
      tNet.begin(), nSamples, nNodes, std::ref(addrWD), std::ref(addrNC), 
      std::ref(addrCV),modsPresent, modNodePresentMap, modIdxMap, nullIdx, 
      nullMap, nulls.memptr(), nPerm, chunkPerms.at(ii), startIdx.at(ii), 
      progress.memptr(), nThreads, ii, std::ref(interrupted)
    );
  }

  MonitorProgress(nPerm, progress.memptr(), nThreads, interrupted, verboseFlag);

  // Wait for all the threads to finish
  for (unsigned int ii = 0; ii < nThreads; ++ii) {
    tt[ii].join();
  }

  // Convert any NaNs or Infinites to NA_REALs
  nulls.elem(arma::find_nonfinite(nulls)).fill(NA_REAL);
  obs.elem(arma::find_nonfinite(obs)).fill(NA_REAL);
  
  // Construct permutation names
  std::vector<std::string> permNames(nPerm);
  for (unsigned int ii = 0; ii < permNames.size(); ++ii) {
    permNames[ii] = "permutation." + std::to_string(ii + 1);
  }
  
  // Convert cube of null distribution objects into an R array before returning
  Rcpp::NumericVector nullsArray (nulls.begin(), nulls.end());
  nullsArray.attr("dim") = Rcpp::IntegerVector::create(mods.size(), 7, nPerm);
  nullsArray.attr("dimnames") = Rcpp::List::create(
    modules, Rcpp::CharacterVector(statnames.begin(), statnames.end()), 
    permNames);
  
  // Convert matrix of observed test statistics into an R object before
  // returning
  Rcpp::NumericMatrix observed (obs.n_rows,  obs.n_cols, obs.begin());
  colnames(observed) = Rcpp::CharacterVector(statnames.begin(), statnames.end());
  rownames(observed) = modules;
  
  return Rcpp::List::create(
    Rcpp::Named("nulls") = nullsArray,
    Rcpp::Named("observed") = observed
  );
}
