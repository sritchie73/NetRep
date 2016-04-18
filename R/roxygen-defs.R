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
#'   backend. If \code{NULL} (default) the all but one core on the machine will
#'   be used.
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
#' @param orderNodesBy \code{NULL} (default), \code{NA}, or a vector of dataset
#'  names or indices. Controls how nodes are ordered on the plot (see details).
#' @param orderSamplesBy \code{NULL} (default), \code{NA}, or a vector 
#'  containing a single dataset name or index. Controls how samples are ordered 
#'  on the plot (see details).
#' @param plotNodeNames logical; controls whether the node names are 
#'  drawed on the bottom axis.
#' @param plotSampleNames logical; controls whether the sample names are 
#'  drawed on the left axis.
#' @param plotModuleNames logical; controls whether module names are drawed.
#'  The default is for module names to be drawed when multiple \code{modules} 
#'  are drawn.
#' @param main title for the plot.
#' @param main.line the number of lines into the top margin at which the plot
#'  title will be drawn.  
#' @param drawBorders logical; if \code{TRUE}, borders are drawn around the 
#'  \emph{weighted degree}, \emph{node conribution}, and \emph{module summary}
#'  bar plots.
#' @param lwd line width for borders and axes.
#' @param naxt.line the number of lines into the bottom margin at which the node
#'  names will be drawn.
#' @param saxt.line the number of lines into the left margin at which the sample
#'  names will be drawn.
#' @param maxt.line the number of lines into the bottom margin at which the 
#'  module names will be drawn.
#' @param xaxt.line the number of lines into the bottom margin at which the 
#'  x-axis tick labels will be drawn on the module summary bar plot.
#' @param xaxt.tck the size of the x-axis ticks for the module summary bar 
#'  plot.
#' @param xlab.line the number of lines into the bottom margin at which the 
#'  x axis label on the \emph{module summary} bar plot(s) will be drawn. 
#' @param yaxt.line the number of lines into the left margin at which the 
#'  y-axis tick labels will be drawn on the weighted degree and node 
#'  contribution bar plots. 
#' @param ylab.line the number of lines into the left margin at which the 
#'  y axis labels on the \emph{weighted degree} and \emph{node contribution} 
#'  bar plots will be drawn. 
#' @param yaxt.tck the size of the y-axis ticks for the weighted degree and 
#'  node contribution bar plots.
#' @param laxt.line the distance from the legend to draw the legend axis 
#'  labels, as multiple of \code{laxt.tck}.
#' @param laxt.tck size of the ticks on each axis legend relative to the
#'  size of the correlation, edge weights, and data matrix heatmaps.
#' @param legend.main.line the distance from the legend to draw the legend 
#'  title.
#' @param cex.axis relative size of the node and sample names.
#' @param cex.lab relative size of the module names and legend titles.
#' @param cex.main relative size of the plot titles.
#' @param dataCols a character vector of colors to create a gradient from for
#'  the data heatmap (see details). Automatically determined if \code{NA} or 
#'  \code{NULL}.
#' @param dataRange the range of values to map to the \code{dataCols} gradient
#'  (see details). Automatically determined if \code{NA} or \code{NULL}.
#' @param corCols a character vector of colors to create a gradient from for
#'  the correlation structure heatmap (see details).
#' @param corRange the range of values to map to the \code{corCols} gradient
#'  (see details).
#' @param netCols a character vector of colors to create a gradient from for
#'  the network edge weight heatmap (see details).
#' @param netRange the range of values to map to the \code{corCols} gradient
#'  (see details). Automatically determined if \code{NA} or \code{NULL}.
#' @param degreeCol color to use for the weighted degree bar plot.
#' @param contribCols color(s) to use for the node contribution bar plot 
#'  (see details).
#' @param summaryCols color(s) to use for the node contribution bar plot 
#'  (see details).
#' @param naCol color to use for missing nodes and samples on the data, 
#'  correlation structure, and network edge weight heat maps.
#' @param dryRun logical; if \code{TRUE}, only the axes and labels will be 
#'  drawed.
#' 
#' @name plot_params
NULL

