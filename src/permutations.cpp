#include <thread>
#include "utils.hpp"
#include "netStats.hpp"

using namespace std;
using namespace arma;
using namespace Rcpp;

//' Calculate the observed test statistics
//'
//' @export
// [[Rcpp::export]]
NumericMatrix PermutationProcedure (
  NumericMatrix dData, NumericMatrix dCorr, NumericMatrix dNet,
  NumericMatrix tData, NumericMatrix tCorr, NumericMatrix tNet,
  CharacterVector moduleAssignments, CharacterVector modules
) {
  /* First, we need to create pointers to the memory holding each
   * NumericMatrix that can be recognised by the armadillo library.
   */
  const mat& dDataPtr = arma::mat(
    dData.begin(), dData.nrow(), dData.ncol(), false, true
  );
  const mat& dCorrPtr = arma::mat(
    dCorr.begin(), dCorr.nrow(), dCorr.ncol(), false, true
  );
  const mat& dNetPtr = arma::mat(
    dNet.begin(), dNet.nrow(), dNet.ncol(), false, true
  );
  const mat& tDataPtr = arma::mat(
    tData.begin(), tData.nrow(), tData.ncol(), false, true
  );
  const mat& tCorrPtr = arma::mat(
    tCorr.begin(), tCorr.nrow(), tCorr.ncol(), false, true
  );
  const mat& tNetPtr = arma::mat(
    tNet.begin(), tNet.nrow(), tNet.ncol(), false, true
  );
  
  /* Next we will scale the matrix data
   */
  const mat dDataScaled = Scale(dDataPtr);
  const mat tDataScaled = Scale(tDataPtr);

  /* Next, we need to create three mappings:
   *  - From node IDs to indices in the discovery dataset
   *  - From node IDs to indices in the test dataset
   *  - From modules to node IDs
   */
  const namemap dIdxMap = makeIdxMap(as<vector<string>>(colnames(dNet)));
  const namemap tIdxMap = makeIdxMap(as<vector<string>>(colnames(tNet)));
  const stringmap modNodeMap = MakeModMap(moduleAssignments);

  // What modules do we actually want to analyse?
  const vector<string> mods = as<vector<string>>(modules);
  const namemap modIdxMap = makeIdxMap(mods);
  
  /**
   * CODE FROM THIS POINT ONWARDS CANNOT USE RCPP!!
   * 
   * The R API is single threaded, we *must not* access it
   * at all during our parallel calls.
   *  
   **/

  // We will save the weighted degree, node contribution, and correlation
  // vectors in the observed dataset so we don't have to compute these at
  // every permutation.
  vecmap obsWD, obsNC, obsCV;
  for (auto mi = mods.begin(); mi != mods.end(); ++mi) {
    // Get the node indices in the discovery dataset for this module
    string mod = *mi;
    uvec dIdx = GetNodeIdx(mod, modNodeMap, dIdxMap);

    // Calculate the network properties and insert into their storage containers
    obsCV.emplace(mod, CorrVector(dCorrPtr, dIdx));

    uvec dRank = sortNodes(dIdx); // Sort nodes for sequential memory access
    obsWD.emplace(mod, WeightedDegree(dNetPtr, dIdx)(dRank));
    vec dSP = SummaryProfile(dDataScaled, dIdx);
    obsNC.emplace(mod, NodeContribution(dDataScaled, dIdx, dSP)(dRank));
  }

  /* Now calculate the observed test statistics
   *
   */
  mat obs (mods.size(), 7);
  for (auto mi = mods.begin(); mi != mods.end(); ++mi) {
    
    // What module are we analysing, and what index does it have internally?
    string mod = *mi;
    unsigned int modIdx = modIdxMap.at(mod);
    
    // Get the node indices in the test dataset for this module
    uvec tIdx = GetNodeIdx(mod, modNodeMap, tIdxMap);

    // Now calculate required properties in the test dataset
    vec tCV = CorrVector(tCorrPtr, tIdx);
    
    uvec tRank = sortNodes(tIdx); // Sort nodes for sequential memory access
    vec tWD = WeightedDegree(tNetPtr, tIdx)(tRank);
    vec tSP = SummaryProfile(tDataScaled, tIdx);
    vec tNC = NodeContribution(tDataScaled, tIdx, tSP)(tRank);

    /* Calculate and store test statistics in the appropriate location in the 
     * results matrix
     */
    obs(modIdx, 0) = AverageEdgeWeight(tWD);
    obs(modIdx, 1) = ModuleCoherence(tNC);
    obs(modIdx, 2) = Correlation(obsCV[mod], tCV);
    obs(modIdx, 3) = Correlation(obsWD[mod], tWD);
    obs(modIdx, 4) = Correlation(obsNC[mod], tNC);
    obs(modIdx, 5) = SignAwareMean(obsCV[mod], tCV);
    obs(modIdx, 6) = SignAwareMean(obsNC[mod], tNC);
  }
  
  // Run permutation procedure using multiple threads
  // int nThreads = nCores[0]
  // 
  // if (nThreads == 1) {
  //   testFun(0, dDataPtr, dCorrPtr, dNetPtr, tDataPtr, tCorrPtr, tNetPtr,
  //           res1, res2, res3, res4, res5, res6);
  // } else {
  //   std::thread *tt = new std::thread[nThreads];
  //   
  //   for (int i = 0; i < nThreads; ++i) {
  //     tt[i] = std::thread(
  //       testFun, i, std::ref(dDataPtr), std::ref(dCorrPtr), std::ref(dNetPtr), 
  //       std::ref(tDataPtr), std::ref(tCorrPtr), std::ref(tNetPtr),
  //       std::ref(res1), std::ref(res2), std::ref(res3),
  //       std::ref(res4), std::ref(res5), std::ref(res6)
  //     );
  //   }
  //   
  //   for (int i = 0; i < nThreads; ++i) {
  //     tt[i].join();
  //   }
  // }
  
  /**
  * RCPP ALLOWED
  **/
  
  /* Convert matrix of observed test statistics into an R object before
   * returning
   */
  NumericMatrix observed (obs.n_rows,  obs.n_cols, obs.begin());
  colnames(observed) = CharacterVector::create("avg.weight", "coherence", 
           "cor.cor", "cor.degree", "cor.contrib", "avg.cor", "avg.contrib");
  rownames(observed) = modules;
  
  return observed;
}
