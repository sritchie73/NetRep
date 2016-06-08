#include "utils.h"
#include "netStats.h"

//' Calculate the network properties 
//' 
//' @details
//' \subsection{Input expectations:}{
//'   Note that this function expects all inputs to be sensible, as checked by
//'   the R function 'checkUserInput' and processed by 'networkProperties'. 
//'   
//'   These requirements are:
//'   \itemize{
//'   \item{The ordering of node names across 'data' and 'net' is consistent.}
//'   \item{The columns of 'data' are the nodes.}
//'   \item{'net' is a square matrix, and its rownames are identical to its 
//'         column names.}
//'   \item{'moduleAssigments' is a named character vector, where the names
//'         represent node labels found in the discovery dataset. Unlike 
//'         'PermutationProcedure', these may include nodes that are not 
//'         present in 'data' and 'net'.}
//'   \item{The module labels specified in 'modules' must occur in 
//'         'moduleAssignments'.}
//'   }
//' }
//' 
//' @param data data matrix from the dataset in which to calculate the network
//'   properties.
//' @param net adjacency matrix of network edge weights between all pairs of 
//'   nodes in the dataset in which to calculate the network properties.
//' @param moduleAssignments a named character vector containing the module 
//'   each node belongs to in the discovery dataset. 
//' @param modules a character vector of modules for which to calculate the 
//'   network properties for.
//' 
//' @return a list containing the summary profile, node contribution, module
//'   coherence, weighted degree, and average edge weight for each 'module'.
// [[Rcpp::export]]
Rcpp::List NetworkProperties (
    Rcpp::NumericMatrix data, Rcpp::NumericMatrix net, 
    Rcpp::CharacterVector moduleAssignments,
    Rcpp::CharacterVector modules
) {
  // First, scale the matrix data
  unsigned int nSamples = data.nrow();
  unsigned int nNodes = data.ncol();
  arma::mat scaledData = Scale(data.begin(), nSamples, nNodes);
  
  R_CheckUserInterrupt(); 
  
  // convert the colnames / rownames to C++ equivalents
  const std::vector<std::string> nodeNames (Rcpp::as<std::vector<std::string>>(colnames(net)));
  const std::vector<std::string> sampleNames (Rcpp::as<std::vector<std::string>>(rownames(data)));
  
  /* Next, we need to create two mappings:
  *  - From node IDs to indices in the dataset of interest
  *  - From modules to node IDs
  *  - From modules to only node IDs present in the dataset of interest
  */
  const namemap nodeIdxMap = MakeIdxMap(nodeNames);
  const stringmap modNodeMap = MakeModMap(moduleAssignments);
  const stringmap modNodePresentMap = MakeModMap(moduleAssignments, nodeIdxMap);
  
  // What modules do we actually want to analyse?
  const std::vector<std::string> mods (Rcpp::as<std::vector<std::string>>(modules));
  
  R_CheckUserInterrupt(); 
  
  // Calculate the network properties for each module
  std::string mod; // iterators
  unsigned int mNodesPresent, mNodes;
  arma::uvec nodeIdx, propIdx, nodeRank;
  namemap propIdxMap;
  std::vector<std::string> modNodeNames; 
  arma::vec WD, SP, NC; // results containers
  double avgWeight, coherence; 
  Rcpp::NumericVector degree, summary, contribution; // for casting to R equivalents
  Rcpp::List results; // final storage container
  for (auto mi = mods.begin(); mi != mods.end(); ++mi) {
    // What nodes are in this module?
    mod = *mi;
    modNodeNames = GetModNodeNames(mod, modNodeMap);
    
    // initialise results containers with NA values for nodes not present in
    // the dataset we're calculating the network properties in.
    degree = Rcpp::NumericVector(modNodeNames.size(), NA_REAL);
    contribution = Rcpp::NumericVector(modNodeNames.size(), NA_REAL);
    summary = Rcpp::NumericVector(nSamples, NA_REAL);
    avgWeight = NA_REAL;
    coherence = NA_REAL;
    degree.names() = modNodeNames;
    contribution.names() = modNodeNames;
    
    // Create a mapping between node names and the result vectors
    propIdxMap = MakeIdxMap(modNodeNames);
    
    // Get just the indices of nodes that are present in the requested dataset
    nodeIdx = GetNodeIdx(mod, modNodePresentMap, nodeIdxMap);
    mNodesPresent = nodeIdx.n_elem;
    
    // And a mapping of those nodes to the initialised vectors
    propIdx = GetNodeIdx(mod, modNodePresentMap, propIdxMap);
    mNodes = propIdx.n_elem;
    
    // Calculate the properties if the module has nodes in the test dataset
    if (nodeIdx.n_elem > 0) {
      // sort the node indices for sequential memory access
      nodeRank = SortNodes(nodeIdx.memptr(), mNodesPresent);
      
      WD = WeightedDegree(net.begin(), nNodes, nodeIdx.memptr(), mNodesPresent);
      WD = WD(nodeRank); // reorder results
      
      avgWeight = AverageEdgeWeight(WD.memptr(), WD.n_elem);
      R_CheckUserInterrupt(); 
      
      SP = SummaryProfile(scaledData.memptr(), nSamples, nNodes, 
                          nodeIdx.memptr(), mNodesPresent);
      R_CheckUserInterrupt(); 
      
      NC = NodeContribution(scaledData.memptr(), nSamples, nNodes, 
                            nodeIdx.memptr(), mNodesPresent, SP.memptr());
      NC = NC(nodeRank); // reorder results
      
      coherence = ModuleCoherence(NC.memptr(), mNodesPresent);
      R_CheckUserInterrupt(); 
      
      // Fill the results vectors appropriately
      Fill(degree, WD.memptr(), mNodesPresent, propIdx.memptr(), mNodes);
      Fill(contribution, NC.memptr(), mNodesPresent, propIdx.memptr(), mNodes);
      summary = Rcpp::NumericVector(SP.begin(), SP.end());
    }
    summary.names() = sampleNames;
    
    results.push_back(
      Rcpp::List::create(
        Rcpp::Named("summary") = summary, 
        Rcpp::Named("contribution") = contribution, 
        Rcpp::Named("coherence") = coherence, 
        Rcpp::Named("degree") = degree,
        Rcpp::Named("avgWeight") = avgWeight
      )
    );
  }
  results.names() = mods;

  return(results);
}

//' Calculate the network properties, data matrix not provided
//' 
//' @details
//' \subsection{Input expectations:}{
//'   Note that this function expects all inputs to be sensible, as checked by
//'   the R function 'checkUserInput' and processed by 'networkProperties'. 
//'   
//'   These requirements are:
//'   \itemize{
//'   \item{'net' is a square matrix, and its rownames are identical to its 
//'         column names.}
//'   \item{'moduleAssigments' is a named character vector, where the names
//'         represent node labels found in the discovery dataset. Unlike 
//'         'PermutationProcedure', these may include nodes that are not 
//'         present in 'data' and 'net'.}
//'   \item{The module labels specified in 'modules' must occur in 
//'         'moduleAssignments'.}
//'   }
//' }
//' 
//' @param net adjacency matrix of network edge weights between all pairs of 
//'   nodes in the dataset in which to calculate the network properties.
//' @param moduleAssignments a named character vector containing the module 
//'   each node belongs to in the discovery dataset. 
//' @param modules a character vector of modules for which to calculate the 
//'   network properties for.
//' 
//' @return a list containing the summary profile, node contribution, module
//'   coherence, weighted degree, and average edge weight for each 'module'.
// [[Rcpp::export]]
Rcpp::List NetworkPropertiesNoData (
    Rcpp::NumericMatrix net, 
    Rcpp::CharacterVector moduleAssignments,
    Rcpp::CharacterVector modules
) {
  // convert the colnames / rownames to C++ equivalents
  const std::vector<std::string> nodeNames (Rcpp::as<std::vector<std::string>>(colnames(net)));
  unsigned int nNodes = net.ncol();
  
  R_CheckUserInterrupt(); 
  
  /* Next, we need to create two mappings:
  *  - From node IDs to indices in the dataset of interest
  *  - From modules to node IDs
  *  - From modules to only node IDs present in the dataset of interest
  */
  const namemap nodeIdxMap = MakeIdxMap(nodeNames);
  const stringmap modNodeMap = MakeModMap(moduleAssignments);
  const stringmap modNodePresentMap = MakeModMap(moduleAssignments, nodeIdxMap);
  
  // What modules do we actually want to analyse?
  const std::vector<std::string> mods (Rcpp::as<std::vector<std::string>>(modules));
  
  R_CheckUserInterrupt(); 
  
  // Calculate the network properties for each module
  std::string mod; // iterators
  unsigned int mNodesPresent, mNodes;
  arma::uvec nodeIdx, propIdx, nodeRank;
  namemap propIdxMap;
  std::vector<std::string> modNodeNames; 
  arma::vec WD; // results containers
  double avgWeight; 
  Rcpp::NumericVector degree; // for casting to R equivalents
  Rcpp::List results; // final storage container
  for (auto mi = mods.begin(); mi != mods.end(); ++mi) {
    // What nodes are in this module?
    // modNodeNames = names(moduleAssignments[moduleAssignments == mod])
    mod = *mi;
    modNodeNames = GetModNodeNames(mod, modNodeMap);
    
    // initialise results containers with NA values for nodes not present in
    // the dataset we're calculating the network properties in.
    degree = Rcpp::NumericVector(modNodeNames.size(), NA_REAL);
    avgWeight = NA_REAL;
    degree.names() = modNodeNames;
    
    // Create a mapping between node names and the result vectors
    propIdxMap = MakeIdxMap(modNodeNames);
    
    // Get just the indices of nodes that are present in the requested dataset
    nodeIdx = GetNodeIdx(mod, modNodePresentMap, nodeIdxMap);
    mNodesPresent = nodeIdx.n_elem;
    
    // And a mapping of those nodes to the initialised vectors
    propIdx = GetNodeIdx(mod, modNodePresentMap, propIdxMap);
    mNodes = propIdx.n_elem;

    // Calculate the properties if the module has nodes in the test dataset
    if (nodeIdx.n_elem > 0) {
      // sort the node indices for sequential memory access
      nodeRank = SortNodes(nodeIdx.memptr(), mNodesPresent);
      
      WD = WeightedDegree(net.begin(), nNodes, nodeIdx.memptr(), mNodesPresent);
      WD = WD(nodeRank); // reorder results
      
      avgWeight = AverageEdgeWeight(WD.memptr(), WD.n_elem);
      R_CheckUserInterrupt(); 
      
      // Fill the results vectors appropriately
      Fill(degree, WD.memptr(), mNodesPresent, propIdx.memptr(), mNodes);
    }

    results.push_back(
      Rcpp::List::create(
        Rcpp::Named("degree") = degree,
        Rcpp::Named("avgWeight") = avgWeight
      )
    );
  }
  results.names() = mods;
  
  return(results);
}
