#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param data a list of numeric matrices. Each entry of the list corresponds to
#'   a dataset and contains the data used to infer the interaction network
#'   between variables (e.g. genes). Expects matrix columns to correspond to
#'   variables and matrix rows to correspond to samples.
#' @param correlation a list of matrices. Each entry of the list corresponds to a 
#'  dataset and contains an \eqn{n * n} matrix of the correlation between 
#'  each pair of variables in the dataset.
#' @param network a list of matrices. Each entry of the list corresponds to a 
#'  dataset and contains an \eqn{n * n} matrix of the network edge weights 
#'  between each pair of variables in the dataset.
#' @param moduleAssignments a vector containing the module each variable belongs
#'  to in the discovery dataset. If there are multiple discovery datasets 
#'  then this argument should be a list of such vectors.
#' @param modules a list of vectors, one for each \code{discovery} dataset, 
#'  of modules to perform the analysis on. The default is to analyse all modules
#'  with the exception of those specified in \code{backgroundLabel}.
#' @param backgroundLabel a single label that nodes that do not belong to any
#'  module are assigned. The default is "0".
#' @param discovery a vector of names or indices denoting the discovery dataset(s).
#' @param test a list of vectors of names or indices denoting the test datasets
#'  for each \code{discovery} dataset. Alternatively can be provided as vector
#'  if the test datasets are the same for all 'discovery' datasets (e.g. for 
#'  performing a pairwise comparison).
#' @param verbose logical; should progress be reported? Default is \code{TRUE}.
#' 
#' @name common_params
NULL

#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param simplify logical; if \code{TRUE}, simplify the structure of the output
#'  list if possible (see Return Value).
#'  
#' @name simplify_param
NULL

#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param nCores number of cores to parallelise the calculation of network 
#'  properties over. Ignored if the user has already registered a parallel 
#'  backend.
#'
#' @name par_param
NULL

#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param orderModules logical; if \code{TRUE} modules ordered by clustering 
#'   their summary vectors. If \code{FALSE} modules are returned in the order
#'   provided.
#'
#' @name orderModules_param
NULL

#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param orderNodesBy one of "discovery", "test" or "none". Controls how nodes
#' are ordered on the plot (see details).
#' @param orderSamplesBy one of "discovery", "test" or "none". Controls how 
#' samples are ordered on the plot (see details).
#' @param plotNodeNames logical; controls whether the node names are 
#'  rendered on the bottom axis.
#' @param plotSampleNames logical; controls whether the sample names are 
#'  rendered on the left axis.
#' @param plotModuleNames logical; controls whether module names are rendered.
#'  The default is for module names to be rendered when multiple \code{modules} 
#'  are drawn.
#' @param drawBorders logical; if \code{TRUE}, borders are drawn around the 
#'  connectivity, module membership, and module summary bar plots.
#' @param border.width line width for borders.
#' @param gaxt.line the number of lines into the bottom margin at which the node
#'  names will be drawn.
#' @param saxt.line the number of lines into the left margin at which the sample
#'  names will be drawn.
#' @param maxt.line the number of lines into the bottom margin at which the 
#'  module names will be drawn.
#' @param legend.tick.size size of the ticks on each axis legend relative to the
#'  size of the correlation, edge weights, and data matrix heatmaps.
#' @param laxt.line the distance from the legend to render the legend axis 
#'  labels, as multiple of \code{legend.tick.size}.
#' @param cex.axis relative size of the node and sample names.
#' @param cex.lab relative size of the module names and legend titles.
#' @param cex.main relative size of the plot titles.
#' 
#' @name plot_params
NULL

