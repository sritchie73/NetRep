#' Template parameters
#' 
#' Template parameters to be imported into other function documentation. This 
#' is not intended to be a stand-alone help file.
#'
#' @param data a list of matrices, one for each dataset. Each 
#'  entry of the list should be a numeric \code{\link{bigMatrix}} of data used
#'  to infer the interaction network for that dataset. The columns should
#'  correspond to variables in the data (nodes in the network) and rows to 
#'  samples in that dataset.
#' @param correlation a list of matrices, one for each dataset. Each entry of
#'   the list should be a \eqn{n * n} numeric \code{\link{bigMatrix}} where each
#'   element contains the correlation coefficient between nodes \eqn{i} and
#'   \eqn{j} in data used to infer the interaction network for that dataset.
#' @param network a list of interaction networks, one for each dataset. Each 
#'  entry of the list should be a \eqn{n * n} numeric \code{\link{bigMatrix}}
#'  where each element contains the edge weight between nodes \eqn{i} and 
#'  \eqn{j} in the inferred network for that dataset. 
#' @param moduleAssignments a list of vectors, one for each \emph{discovery} 
#'   dataset, containing the module assignments for each node in that dataset.
#' @param modules a list of vectors, one for each \code{discovery} dataset, 
#'  of modules to perform the analysis on. If unspecified, all modules
#'  in each \code{discovery} dataset will be analysed, with the exception of 
#'  those specified in \code{backgroundLabel} argument.
#' @param backgroundLabel a single label given to nodes that do not belong to 
#'  any module in the \code{moduleAssignments} argument.
#' @param discovery a vector of names or indices denoting the \emph{discovery}
#'  dataset(s) in the \code{data}, \code{correlation}, \code{network}, 
#'  \code{moduleAssignments}, \code{modules}, and \code{test} lists. 
#' @param test a list of vectors, one for each \code{discovery} dataset,
#'  of names or indices denoting the \emph{test} dataset(s) in the \code{data}, 
#'  \code{correlation}, and \code{network} lists.
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
#'   properties over. Ignored if the user has already registered a parallel 
#'   backend.If \code{NULL} (default) the maximum number of cores on the machine
#'   will be used.
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

