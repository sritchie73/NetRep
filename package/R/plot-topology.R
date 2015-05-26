#' Plot a topological feature of network module
#' 
#' Functions for plotting the topology of a network module.
#' 
#' @param geneExpression optional; a list of 
#'   \code{\link[=bigMatrix-class]{bigMatrix}} objects, each containing the gene
#'   expression data for a datset of interest (see details). Columns are
#'   expected to be genes, rows samples. If not provided, the expression,
#'   module membership, and summary expression will not be plotted.
#' @param coexpression a list of 'bigMatrix' objects, each containing the gene
#'   coexpression for a dataset of interest (see details).
#' @param adjacency a list of 'bigMatrix' objects, each containing the gene
#'   adjacencies for a dataset of interest (see details).
#' @param moduleAssignments a list of named vectors assigning genes to modules 
#'   in each dataset of interest (see details).
#' @param modules a vector of modules to apply the function to (see details).
#' @param discovery name or index denoting which dataset the module of interest
#'   was discovered in (see details).
#' @param test name or index denoting which dataset to apply the function to 
#'   (see details).
#' @param symmetric logical; controls whether the coexpression and adjacency 
#'  heatmaps are drawn as symmetric (square) heatmaps or asymettric triangle 
#'  heatmaps. If symmetric, then the gene and module names will also be rendered
#'  on the left axis.
#' @param orderGenesBy one of "discovery", "test", or "none" indicating which
#'  dataset the function should use to order the genes (see details).
#' @param orderSamplesBy one of "discovery", "test", or "none" indicating which
#'  dataset the function should use to order the samples (see details).
#' @param orderModules logical; if \code{TRUE} modules will be ordered by the 
#'  similarity of their summary expression profiles (see details). The default
#'  is \code{TRUE} if gene expression for the \code{test} dataset is provided.
#' @param plotGeneNames logical; controls whether the gene names are 
#'  rendered on the bottom axis when using \code{plotExpression},
#'  \code{plotCoxpression}, \code{plotAdjacency}, \code{plotConnectivity},
#'  or \code{plotModuleMembership}.
#' @param plotSampleNames logical; controls whether the sample names are 
#'  rendered on the left axis when using \code{plotSummaryExpression} or 
#'  \code{plotExpression}.
#' @param plotModuleNames logical; controls whether module names are rendered 
#'  above the bar charts when using \code{plotSummaryExpression}, or on the 
#'  bottom axis when using \code{plotExpression}, \code{plotCoexpression},
#'  \code{plotAdjacency}, \code{plotConnectivity}, or 
#'  \code{plotModuleMembership}. By default, module names are only rendered 
#'  when multiple \code{modules} are specified.
#' @param palette a vector of colors to use for each plot (see details).
#' @param drawBorders logical; if \code{TRUE}, borders are drawn around the bars
#'  in \code{plotModuleMembership}, \code{plotConnectivity}, and
#'  \code{plotSummaryExpression}.
#' @param plotLegend logical; controls whether a legend is drawn when using
#'  \code{plotCoexpression}, \code{plotAdjacency}, or \code{plotExpression}.
#' @param gaxt.line the number of lines into the bottom margin at which the gene
#'  names will be drawn.
#' @param saxt.line the number of lines into the left margin at which the sample
#'  names will be drawn.
#' @param maxt.line the number of lines into the bottom margin at which the 
#'  module names will be drawn.
#' @param legend.main title for the legend.
#' @param main title for each plot.
#' @param legend.position the distance from the plot to start the legend, as a
#'  proportion of the plot width.
#' @param legend.tick.size size of the ticks on the axis legend.
#' @param laxt.line the distance from the legend to render the legend axis 
#'  labels, as multiple of \code{legend.tick.size}.
#' @param cex.axis relative size of the gene and sample names.
#' @param cex.lab relative size of the module names and legend titles.
#' @param cex.main relative size of the plot titles.
#' @param horizontal logical; controls whether the legend is rendered 
#'  horizontally or vertically when using \code{plotExpressionLegend},
#'  \code{plotCoexpressionLegend} or \code{plotAdjacencyLegend}.
#'
#' @details
#'  \subsection{Input data structure:}{
#'   This function allows for input data formatted in a number of ways. Where 
#'   there are multiple datasets of interest (e.g. multiple tissues, or a 
#'   discovery dataset and an independent test dataset) the arguments 
#'   \code{geneExpression}, \code{coexpression}, and \code{adjacency} should be
#'   \code{\link[=list]{lists}} where each element contains the matrix data for 
#'   each respective dataset. This matrix data should be stored as a 'bigMatrix'
#'   object (see \link[=bigMatrix-get]{converting matrix data to 'bigMatrix'
#'   data}). Alternatively, if only one dataset is of interest, the
#'   \code{geneExpression}, \code{coexpression}, and \code{adjacency} arguments
#'   will also each accept a single 'bigMatrix' object.
#'   
#'   Similarly, the \code{moduleAssignments} argument expects a list of named
#'   vectors, which contain the the module assignments for each gene in the
#'   respective dataset. List elements corresponding to datasets where module
#'   discovery has not been performed should contain \code{NULL}, unless the
#'   datasets are named throughout the function arguments. I.e. where the
#'   \code{\link{names}} of \code{geneExpression}, \code{coexpression}, and
#'   \code{adjacency} correspond to the names of each dataset of interest, the
#'   names of the \code{discovery} dataset can be used to look up the respective
#'   module assignments in the \code{moduleAssignments} list. If module
#'   discovery has only been performed in one dataset, then the
#'   \code{moduleAssignments} will also accept a named vector.
#'   
#'   The \code{discovery} arguments specifies which dataset the \code{modules} 
#'   of interest were discovered in, and the \code{test} argument specifies 
#'   which dataset to plot those module(s) in. The \code{orderGenesBy} and
#'   \code{orderSamplesBy} arguments control how genes and samples are ordered
#'   on the plot. These arguments are ignored if data is provided for only one 
#'   dataset.
#' }
#' \subsection{'bigMatrix' vs. 'matrix' input data:}{
#'   Although the function expects \code{\link[=bigMatrix-class]{bigMatrix}}
#'   data, regular 'matrix' objects are also accepted. In this case, the
#'   'matrix' data is temporarily converted to 'bigMatrix' by the function. This
#'   conversion process involves writing out each matrix as a binary file on
#'   disk, which can take a long time for large datasets. It is strongly
#'   recommended for the user to store their data as 'bigMatrix' objects, as the
#'   \link{modulePreservation} function, \link{networkProperties} function,
#'   \link[=plotModule]{plotting} \link[=plotTopology]{functions},
#'   \link[=geneOrder]{gene} and \link[=sampleOrder]{sample} ordering also
#'   expect 'bigMatrix' objects. Further, 'bigMatrix' objects have a number of
#'   benefits, including instantaneous load time from any future R session, and
#'   parallel access from mutliple independent R sessions. Methods are provided
#'   for \link[=bigMatrix-get]{converting to, loading in}, and 
#'   \link[=bigMatrix-out]{writing out} 'bigMatrix' objects.
#' }
#' \subsection{Gene and sample ordering:}{
#'   By default, genes are ordered in decreasing order of intramodular 
#'   connectivity in the \code{discovery} dataset (see \code{\link{geneOrder}}). 
#'   This facilitates the visual comparison of modules across datasets, as the 
#'   gene ordering will be preserved. Missing genes are colored in grey. This 
#'   behaviour can be change by setting \code{orderGenesBy} to "test", in which
#'   cases genes will be ordered in decreasing order of intramodular 
#'   connectivity in the discovery dataset. Alternatively \code{orderGenesBy}
#'   can be set to "none", in which case genes are rendered in the order they
#'   appear in the discovery dataset.
#'   
#'   When multiple modules are specified, modules are ordered by the similarity
#'   of their summary expression profiles in the \code{orderGenesBy} dataset. 
#'   To disable this behaviour, set \code{orderModules} to \code{FALSE}.
#'   
#'   Sample ordering only applies to \code{plotExpression} and 
#'   \code{plotSummaryExpression}. By default, samples are ordered in descending
#'   order of the summary expression profile for the left-most module appearing
#'   on the plot (see \code{\link{sampleOrder}}. By default, the summary
#'   expression profile is calculated in the \code{test} dataset. This behaviour
#'   can be changed through the \code{orderSamplesBy} argument, however setting 
#'   \code{orderSamplesBy} to "discovery" will only work if samples are present 
#'   in both datasets.
#' }
#' \subsection{Normalised intramodular connectivity:}{
#'   The gene connectivity is normalised by the maximum connectivity in any 
#'   given module when rendered on the bar plot. This facilitates visual 
#'   comparison on genes within a module when multiple modules of differing 
#'   sizes or densities are rendered. Further, although the relative
#'   intramodular connectivity provides information about a genes biological
#'   importance to a module \emph{(1)}, the numeric value is meaningless.
#'   Normalising the connectivity is therefore useful for visual inspection.
#' }
#' \subsection{Customising plot layout:}{
#'   Although reasonable default values for most parameters have been provided,
#'   the rendering of axes and titles may need adjusting depending on the size
#'   of the plot window. The parameters \code{gaxt.line}, \code{saxt.line}, 
#'   \code{maxt.line}, and \code{laxt.line} control the distance from each plot
#'   window that the gene labels, sample labels, module labels, and legend 
#'   labels are rendered. 
#'   
#'   \code{legend.tick.size} controls the length of the 
#'   axis ticks on each of the legends relative to the coexpression, adjacency,
#'   and gene expression plot windows. 
#'   
#'   \code{legend.position} controls the horizontal offset of the legend 
#'   relative to the plot. For the triangle heatmaps, (\code{symmetric=FALSE} in
#'   \code{plotCoexpression} and \code{plotAdjacency}) this controls how far 
#'   left of the plot the legend starts as a proportion of the plot width. For
#'   the square heatmaps (\code{plotExpression}, and \code{symmetric=TRUE} in
#'   \code{plotCoexpression} and \code{plotAdjacency}) this controls how far 
#'   right of the plot the legend starts as a proportion of the plot width.
#'   
#'   \code{cex.main} controls the relative text size of the plot title
#'   (specified by the \code{main} argument). \code{cex.axis} controls the
#'   relative text size of the gene and sample labels. \code{cex.lab} controls
#'   the relative text size of the bar plot axis labels, module labels, and the
#'   legend titles.
#'   
#'   The rendering of gene, sample, and module names can be disabled by setting
#'   \code{plotGeneNames}, \code{plotSampleNames}, and \code{plotModuleNames} to
#'   \code{FALSE}, and the rendering of the legend can be disabled by setting
#'   \code{plotLegend} to \code{FALSE}
#'   
#'   The \code{drawBorders} argument controls whether borders are drawn around
#'   the bars in \code{plotConnectivity}, \code{plotModuleMembership}, and 
#'   \code{plotSummaryExpression}.
#' }
#' \subsection{Customising the color palette:}{
#'   \code{plotCoexpression} and \code{plotCoexpressionLegend} expect the
#'   \code{palette} argument to be a vector of colors to interpolate over when
#'   plotting the coexpression. They expect the first element of the 
#'   \code{palette} vector to be the color used for coexpression values of -1,
#'   and the last element of the \code{palette} vector to be the color used for 
#'   coexpression values of 1.
#'   
#'   \code{plotAdjacency} and \code{plotAdjacencyLegend} expect the
#'   \code{palette} argument to be a vector of colors to interpolate over when
#'   plotting the adjacencies. They expect the first element of the 
#'   \code{palette} vector to be the color used for gene adjacency values of 0,
#'   and the last element of the \code{palette} vector to be the color used for 
#'   coexpression values of 1.
#'   
#'   \code{plotConnectivity} expects \code{palette} to be a single color, a 
#'   vector of colors, one for each gene, or a vector of colors to be repeated.
#'   
#'   \code{plotModuleMembership} expects \code{palette} to be a vector 
#'   containing two colors, the first to be used for genes with negative module
#'   membership values, and the second to be used for genes with positive module
#'   membership values. 
#'   
#'   \code{plotExpression} and \code{plotExpressionLegend} expect the 
#'   \code{palette} argument to be a vector of colors to interpolate over when 
#'   plotting the gene expression. In order to accomodate gene expression values
#'   with dramatically different ranges, these functions expect the palette to 
#'   be a diverging set of colors with a centre value corresponding the a gene 
#'   expression value of 0 (e.g. white). The colors used will be balanced around
#'   0: i.e. positive and negative values of gene expression will have the same 
#'   color intensity on a diverging color palette regardless of the actual 
#'   range. For gene expression data not centred around 0 the tail of the vector
#'   will be used: the minimum gene expression value will receive the centre 
#'   palette color, and the last element of the \code{palette} vector will be
#'   used for maximum gene expression value.
#'   
#'   \code{plotSummaryExpression} expects \code{palette} to be a vector 
#'   containing two colors, the first to be used for genes with a negative 
#'   summary expression profile, and the second to be used for genes with a 
#'   positive summary expression profile, regardless of whether the gene 
#'   expression is centred around 0.
#' }
#' 
#' @references
#' \enumerate{
#'    \item{
#'      Langfelder, P., Mischel, P. S. & Horvath, S. \emph{When is hub gene 
#'      selection better than standard meta-analysis?} PLoS One \strong{8}, 
#'      e61505 (2013).
#'    }
#' }
#' 
#' @seealso
#' \code{\link{plotModule}} for a combined plot showing all topological 
#' properties for a network module.
#' 
#' @examples
#' \dontrun{
#' ## Create some example data
#' geA <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geA) <- paste0("Gene_", 1:100)
#' rownames(geA) <- paste0("CohortA_", 1:50)
#' coexpA <- cor(geA) # coexpression
#' adjA <- abs(coexpA)^5 # adjacency
#' moduleAssignments <- sample(1:7, size=100, replace=TRUE)
#' names(moduleAssignments) <- paste0("Gene_", 1:100)
#' 
#' # Create bigMatrix objects for each matrix.
#' geA <- as.bigMatrix(geA, "geA_bm")
#' coexpA <- as.bigMatrix(coexpA, "coexpA_bm")
#' adjA <- as.bigMatrix(adjA, "adjA_bm")
#' 
#' ## Example 1: Plot Module 2 in cohort A.
#' plotExpression(geA, coexpA, adjA, moduleAssignments, modules="2")
#' plotCoexpression(geA, coexpA, adjA, moduleAssignments, modules="2")
#' # alternatively as a square heatmap
#' plotCoexpression(
#'  geA, coexpA, adjA, moduleAssignments, modules="2", symmetric=TRUE
#' )
#' plotAdjacency(geA, coexpA, adjA, moduleAssignments, modules="2")
#' plotConnectivity(geA, coexpA, adjA, moduleAssignments, modules="2")
#' plotModuleMembership(geA, coexpA, adjA, moduleAssignments, modules="2")
#' plotSummaryExpression(geA, coexpA, adjA, moduleAssignments, modules="2")
#' 
#' ## Example 2: Plot an arbitrary set of genes in cohort A
#' plotExpression(geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10])
#' plotCoexpression(geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10])
#' # alternatively as a square heatmap
#' plotCoexpression(
#'  geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10], symmetric=TRUE
#' )
#' plotAdjacency(geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10])
#' plotConnectivity(geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10])
#' plotModuleMembership(geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10])
#' plotSummaryExpression(geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10])
#' 
#' ## Example 3: Plot the topology of two adipose tissue modules in the liver
#' ## tissue data 
#' 
#' geAdipose <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geAdipose) <- paste0("Gene_", 1:100)
#' rownames(geAdipose) <- paste0("Sample_", 1:50)
#' coexpAdipose <- cor(geAdipose) # coexpression
#' adjAdipose <- abs(coexpAdipose)^5 # adjacency
#' adiposeModules <- sample(0:7, size=100, replace=TRUE)
#' names(adiposeModules) <- paste0("Gene_", 1:100)
#' 
#' geLiver <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geLiver) <- paste0("Gene_", 1:100)
#' rownames(geLiver) <- paste0("Sample_", 1:50)
#' coexpLiver <- cor(geLiver) # coexpression
#' adjLiver <- abs(coexpLiver)^6 # adjacency
#' liverModules <- sample(0:12, size=100, replace=TRUE)
#' names(liverModules) <- paste0("Gene_", 1:100)
#'
#' geHeart <- matrix(rnorm(50*100), ncol=100) # gene expression
#' colnames(geHeart) <- paste0("Gene_", 1:100)
#' rownames(geHeart) <- paste0("Sample_", 1:50)
#' coexpHeart <- cor(geHeart) # coexpression
#' adjHeart <- abs(coexpHeart)^4 # adjacency
#' heartModules <- sample(0:5, size=100, replace=TRUE)
#' names(heartModules) <- paste0("Gene_", 1:100)
#' 
#' # Store each input type as a list, where each element corresponds
#' # to one of the datasets
#' geneExpression <- list(
#'   adipose=as.bigMatrix(geAdipose, "geAdipose_bm"),
#'   liver=as.bigMatrix(geLiver, "geLiver_bm"),  
#'   heart=as.bigMatrix(geHeart, "geHeart_bm") 
#' )
#' coexpression <- list(
#'   adipose=as.bigMatrix(coexpAdipose, "coexpAdipose_bm"),
#'   liver=as.bigMatrix(coexpLiver, "coexpLiver_bm"),  
#'   heart=as.bigMatrix(coexpHeart, "coexpHeart_bm") 
#' )
#' adjacency <- list(
#'   adipose=as.bigMatrix(adjAdipose, "adjAdipose_bm"),
#'   liver=as.bigMatrix(adjLiver, "adjLiver_bm"),  
#'   heart=as.bigMatrix(adjHeart, "adjHeart_bm") 
#' )
#' moduleAssignments <- list(
#'   adipose=adiposeModules, liver=liverModules, heart=heartModules
#' )
#' 
#' # Show the plot
#' plotExpression(
#'  geneExpression, coexpression, adjacency, moduleAssignments,
#'  modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' plotCoexpression(
#'  geneExpression, coexpression, adjacency, moduleAssignments,
#'  modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' # alternatively as a square heatmap
#' plotCoexpression(
#'  geneExpression, coexpression, adjacency, moduleAssignments, 
#'  modules=c("3", "7"), discovery="adipose", test="liver", symmetric=TRUE
#' )
#' plotAdjacency(
#'  geneExpression, coexpression, adjacency, moduleAssignments,
#'  modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' plotConnectivity(
#'  geneExpression, coexpression, adjacency, moduleAssignments, 
#'  modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' plotModuleMembership(
#'  geneExpression, coexpression, adjacency, moduleAssignments,
#'  modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' plotSummaryExpression(
#'  geneExpression, coexpression, adjacency, moduleAssignments,
#'  modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' 
#' # clean up bigMatrix files from examples
#' unlink("*_bm*")
#' }
#' 
#' @name plotTopology
NULL

#' \code{plotExpression:} Plot a heatmap of the gene expression for one or more
#' network modules in their discovery dataset, or an independent test dataset.
#' 
#' @rdname plotTopology
#' @export
plotExpression <- function(
  geneExpression, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, orderSamplesBy="test", orderGenesBy="discovery",
  orderModules, plotGeneNames=TRUE, plotSampleNames=TRUE, plotModuleNames,
  main="", palette=expression.palette(), plotLegend=TRUE, 
  legend.main="Expression", gaxt.line=-0.5, saxt.line=-0.5, maxt.line=3, 
  legend.position=0.15, legend.tick.size=0.03, laxt.line=3, cex.axis=0.8, 
  cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  if (is.null(geneExpression))
    stop("Cannot plot gene expression without gene expression data")
  
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderGenesBy <- orderByArgs[pmatch(orderGenesBy, orderByArgs, nomatch=3)]
  orderSamplesBy <- orderByArgs[pmatch(orderSamplesBy, orderByArgs, nomatch=3)]
  
  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- paste0(".temp-objects", getUUID())
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  geneExpression <- unifyDS(dynamicMatLoad(geneExpression, backingpath=tmp.dir))
  coexpression <- unifyDS(dynamicMatLoad(coexpression, backingpath=tmp.dir))
  adjacency <- unifyDS(dynamicMatLoad(adjacency, backingpath=tmp.dir))
  
  # Format optional input data so it doesn't cause cascading error crashes 
  geneExpression <- formatGeneExpression(
    geneExpression, length(coexpression), names(coexpression)
  )
  
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  # Sanity check input for consistency.
  checkSets(
    geneExpression, coexpression, adjacency, moduleAssignments, discovery, test
  )
  
  if (is.null(geneExpression[[test]]))
    stop("Cannot plot gene expression without gene expression data")
  
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  #-----------------------------------------------------------------------------
  # Get ordering of genes and samples in the 'test' dataset by the dataset 
  # specified in 'orderGenesBy' and 'orderSamplesBy'.
  #-----------------------------------------------------------------------------
  props <- networkProperties(
    geneExpression, coexpression, adjacency, moduleAssignments, modules, 
    discovery, test, FALSE
  )
  
  # Determine gene ordering, then sample ordering.
  if (orderGenesBy == "discovery" && discovery != test) {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[discovery]]), FALSE, TRUE)
    
    propsDisc <- networkProperties(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery, FALSE
    ) 
    
    moduleOrder <- names(propsDisc)
    if (length(propsDisc) > 1 && orderModules) {
      # Create a matrix of summary expression profiles to measure the similarity
      seps <- matrix(
        0, ncol=length(propsDisc), 
        nrow=length(propsDisc[[1]]$summaryExpression)
      )
      colnames(seps) <- names(propsDisc)
      for (mi in seq_along(propsDisc)) {
        seps[,mi] <- propsDisc[[mi]]$summaryExpression
      }
      moduleOrder <- names(propsDisc)[hclust(as.dist(1-cor(seps)))$order]
    }
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(sort(
        propsDisc[[mi]]$connectivity, decreasing=TRUE, na.last=TRUE
      ))
    }
  } else if (orderGenesBy == "none") {
    moduleOrder <- names(props)
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(props[[mi]]$connectivity)
    }
  } else {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[test]]), FALSE, TRUE)
    
    # Order modules and samples by the test network
    moduleOrder <- names(props)
    if (length(props) > 1 && orderModules) {
      seps <- matrix(
        0, ncol=length(props), 
        nrow=length(props[[1]]$summaryExpression)
      )
      colnames(seps) <- names(props)
      for (mi in seq_along(props)) {
        seps[,mi] <- props[[mi]]$summaryExpression
      }
      moduleOrder <- names(props)[hclust(as.dist(1-cor(seps)))$order]
    }
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(sort(
        props[[mi]]$connectivity, decreasing=TRUE, na.last=TRUE
      ))
    }
  }
  
  # If test == discovery, we order samples and modules by discovery, and plot
  # the discovery. So we need to make sure the gene expression exists in the
  # discovery
  if (orderSamplesBy == "discovery" && discovery != test) {
    if (is.null(geneExpression[[discovery]])) {
      stop(
        "Expecting gene expression data for the discovery dataset in order",
        " to sort samples"
      )
    }
    if (!exists("propsDisc")) {
      propsDisc <- networkProperties(
        geneExpression, coexpression, adjacency, moduleAssignments, modules,
        discovery, discovery, FALSE
      )
    }
    sampleOrder <- names(sort(
      propsDisc[[1]]$summaryExpression, decreasing=TRUE
    ))
  } else if (orderSamplesBy == "none") {
    sampleOrder <- rownames(geneExpression[[test]])
  } else {
    sampleOrder <- names(sort(
      props[moduleOrder][[1]]$summaryExpression, decreasing=TRUE
    ))
  }
  
  #-----------------------------------------------------------------------------
  # Identify genes and samples from the 'discovery' dataset not present in the 
  # 'test' dataset.
  #-----------------------------------------------------------------------------
  if (all(sampleOrder %nin% rownames(geneExpression[[test]]))) {
    stop(
      "No samples from the 'orderSamplesBy' dataset are present in the",
      " 'test' dataset"
    )
  }
  na.pos.y <- which(sampleOrder %nin% rownames(geneExpression[[test]]))
  if (length(na.pos.y) > 0) {
    presentSamples <- sampleOrder[-na.pos.y]
  } else {
    presentSamples <- sampleOrder
  }
  
  # Handle genes not present in the test dataset
  na.pos.x <- which(geneOrder %nin% colnames(coexpression[[test]]))
  if (length(na.pos.x) > 0) {
    presentGenes <- geneOrder[-na.pos.x]
  } else {
    presentGenes <- geneOrder
  }
  
  #-----------------------------------------------------------------------------
  # Plot the gene expression 
  #-----------------------------------------------------------------------------
  # First we need to set up the color palette for the gene expression, which 
  # includes the fact that the range may be unbalanced around 0.
  ge <- geneExpression[[test]][presentSamples, presentGenes]
  range.ge <- range(ge)
  if (all(range.ge > 0)) {
    palette <- tail(palette, length(palette)/2)
  } else if (all(range.ge < 0)) {
    palette <- head(palette, length(palette)/2)
  } else {
    range.pal <- c(-max(abs(range.ge)), max(abs(range.ge)))
  }
  xaxt <- NULL
  if (plotGeneNames)
    xaxt <- geneOrder
  yaxt <- NULL
  if (plotSampleNames)
    yaxt <- sampleOrder
  plotSquareHeatmap(
    ge, palette, vlim=range.pal, legend.lim=range.ge,
    moduleAssignments[[discovery]][geneOrder], na.pos.x, na.pos.y, 
    xaxt=xaxt, yaxt=yaxt, plotLegend=plotLegend, main=main,
    legend.main=legend.main, plotModuleNames=plotModuleNames,
    xaxt.line=gaxt.line, yaxt.line=saxt.line, legend.tick.size=legend.tick.size,
    laxt.line=laxt.line, legend.line=legend.position, maxt.line=maxt.line
  )
}

#' \code{plotCoexpression:} Plot a heatmap of the gene coexpression for one or more
#' network modules in their discovery dataset, or an independent test dataset.
#' 
#' @rdname plotTopology
#' @export
plotCoexpression <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, symmetric=FALSE, orderGenesBy="discovery", orderModules,
  plotGeneNames=TRUE, plotModuleNames, main="", 
  palette=coexpression.palette(), plotLegend=TRUE, legend.main="Coexpression",
  gaxt.line=-0.5, maxt.line=3, legend.position, legend.tick.size=0.03, 
  laxt.line=2.5, cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  if (class(main) != "character")
    stop("'main' must be a characer vector")

  orderByArgs <- c("discovery", "test", "none")
  orderGenesBy <- orderByArgs[pmatch(orderGenesBy, orderByArgs, nomatch=3)]

  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- paste0(".temp-objects", getUUID())
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  geneExpression <- unifyDS(dynamicMatLoad(geneExpression, backingpath=tmp.dir))
  coexpression <- unifyDS(dynamicMatLoad(coexpression, backingpath=tmp.dir))
  adjacency <- unifyDS(dynamicMatLoad(adjacency, backingpath=tmp.dir))
  
  # Format optional input data so it doesn't cause cascading error crashes 
  geneExpression <- formatGeneExpression(
    geneExpression, length(coexpression), names(coexpression)
  )
  
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  
  #-----------------------------------------------------------------------------
  # Get ordering of genes in the 'test' dataset by the dataset specified in 
  # 'orderGenesBy'.
  #-----------------------------------------------------------------------------
  if (orderGenesBy == "none") {
    geneOrder <- getGenes(moduleAssignments, modules, discovery)
  } else if (orderGenesBy == "discovery") {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[discovery]]), FALSE, TRUE)
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery, FALSE, orderModules
    )
  } else {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[test]]), FALSE, TRUE)
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, test, FALSE, orderModules
    )
  }
 
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  #-----------------------------------------------------------------------------
  # Identify genes from the 'discovery' dataset not present in the 'test' 
  # dataset.
  #-----------------------------------------------------------------------------
  na.pos <- which(geneOrder %nin% colnames(coexpression[[test]]))
  if (length(na.pos) > 0) {
    presentGenes <- geneOrder[-na.pos]
  } else {
    presentGenes <- geneOrder
  }
  
  #-----------------------------------------------------------------------------
  # Plot the gene coexpression 
  #-----------------------------------------------------------------------------
  gaxt <- NULL
  if (plotGeneNames)
    gaxt <- geneOrder
  if (symmetric) {
    if (missing(legend.position))
      legend.position <- 0.2
    plotSquareHeatmap(
      coexpression[[test]][presentGenes, presentGenes], palette, c(-1, 1), 
      moduleAssignments[[discovery]][geneOrder], na.pos, na.pos, 
      xaxt=gaxt, yaxt=gaxt, plotLegend=plotLegend, main=main,
      legend.main=legend.main, plotModuleNames=plotModuleNames,
      xaxt.line=gaxt.line, yaxt.line=gaxt.line, 
      legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
      legend.line=legend.position, maxt.line=maxt.line
    )
  } else {
    if (missing(legend.position))
      legend.position <- 0.1
    plotTriangleHeatmap(
      coexpression[[test]][presentGenes , presentGenes], palette, c(-1, 1),
      moduleAssignments[[discovery]][geneOrder], na.pos, xaxt=gaxt, 
      plotLegend=plotLegend, main=main, legend.main=legend.main, 
      plotModuleNames=plotModuleNames, xaxt.line=gaxt.line,
      legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
      legend.line=legend.position, maxt.line=maxt.line
    )
  }
}

#' \code{plotAdjacency:} Plot a heatmap of the gene adjacencies for one or more
#' network modules in their discovery dataset, or an independent test dataset.
#' 
#' @rdname plotTopology
#' @export
plotAdjacency <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, symmetric=FALSE, orderGenesBy="discovery", orderModules,
  plotGeneNames=TRUE, plotModuleNames, main="", 
  palette=adjacency.palette(), plotLegend=TRUE, legend.main="Adjacency",
  gaxt.line=-0.5, maxt.line=3, legend.position, legend.tick.size=0.03, 
  laxt.line=2.5, cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderGenesBy <- orderByArgs[pmatch(orderGenesBy, orderByArgs, nomatch=3)]
  
  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- paste0(".temp-objects", getUUID())
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  geneExpression <- unifyDS(dynamicMatLoad(geneExpression, backingpath=tmp.dir))
  coexpression <- unifyDS(dynamicMatLoad(coexpression, backingpath=tmp.dir))
  adjacency <- unifyDS(dynamicMatLoad(adjacency, backingpath=tmp.dir))
  
  # Format optional input data so it doesn't cause cascading error crashes 
  geneExpression <- formatGeneExpression(
    geneExpression, length(coexpression), names(coexpression)
  )
  
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  #-----------------------------------------------------------------------------
  # Get ordering of genes in the 'test' dataset by the dataset specified in 
  # 'orderGenesBy'.
  #-----------------------------------------------------------------------------
  if (orderGenesBy == "none") {
    geneOrder <- getGenes(moduleAssignments, modules, discovery)
  } else if (orderGenesBy == "discovery") {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[discovery]]), FALSE, TRUE)
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery, FALSE, orderModules
    )
  } else {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[test]]), FALSE, TRUE)
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, test, FALSE, orderModules
    )
  }
  
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  #-----------------------------------------------------------------------------
  # Identify genes from the 'discovery' dataset not present in the 'test' 
  # dataset.
  #-----------------------------------------------------------------------------
  na.pos <- which(geneOrder %nin% colnames(coexpression[[test]]))
  if (length(na.pos) > 0) {
    presentGenes <- geneOrder[-na.pos]
  } else {
    presentGenes <- geneOrder
  }
  
  #-----------------------------------------------------------------------------
  # Plot the gene adjacency
  #-----------------------------------------------------------------------------
  gaxt <- NULL
  if (plotGeneNames)
    gaxt <- geneOrder
  if (symmetric) {
    if (missing(legend.position))
      legend.position <- 0.2
    plotSquareHeatmap(
      adjacency[[test]][presentGenes, presentGenes], palette, c(0, 1), 
      moduleAssignments[[discovery]][geneOrder], na.pos, na.pos, 
      xaxt=gaxt, yaxt=gaxt, plotLegend=plotLegend, main=main,
      legend.main=legend.main, plotModuleNames=plotModuleNames,
      xaxt.line=gaxt.line, yaxt.line=gaxt.line, 
      legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
      legend.line=legend.position, maxt.line=maxt.line
    )
  } else {
    if (missing(legend.position))
      legend.position <- 0.1
    plotTriangleHeatmap(
      adjacency[[test]][presentGenes , presentGenes], palette, c(0, 1),
      moduleAssignments[[discovery]][geneOrder], na.pos, xaxt=gaxt, 
      plotLegend=plotLegend, main=main, legend.main=legend.main, 
      plotModuleNames=plotModuleNames, xaxt.line=gaxt.line,
      legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
      legend.line=legend.position, maxt.line=maxt.line
    )
  }
}

#' \code{plotModuleMembership:} Plot a bar chart of the module membership for
#' one or more network modules in their discovery dataset, or an independent
#' test dataset.
#' 
#' @rdname plotTopology
#' @export
plotModuleMembership <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, orderGenesBy="discovery", orderModules,
  plotGeneNames=TRUE, plotModuleNames, main="", 
  palette=c("#313695", "#a50026"), drawBorders=FALSE, gaxt.line=-0.5, 
  maxt.line=3, cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  if (is.null(geneExpression))
    stop("Cannot plot module membership without gene expression data")
  
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderGenesBy <- orderByArgs[pmatch(orderGenesBy, orderByArgs, nomatch=3)]
  
  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- paste0(".temp-objects", getUUID())
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  geneExpression <- unifyDS(dynamicMatLoad(geneExpression, backingpath=tmp.dir))
  coexpression <- unifyDS(dynamicMatLoad(coexpression, backingpath=tmp.dir))
  adjacency <- unifyDS(dynamicMatLoad(adjacency, backingpath=tmp.dir))
  
  # Format optional input data so it doesn't cause cascading error crashes
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  if (is.null(geneExpression[[test]]))
    stop("Cannot plot module membership without gene expression data")
  
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  #-----------------------------------------------------------------------------
  # Get ordering of genes in the 'test' dataset by the dataset specified in 
  # 'orderGenesBy'.
  #-----------------------------------------------------------------------------
  # Get the module membership for each module in the test network.
  props <- networkProperties(
    geneExpression, coexpression, adjacency, moduleAssignments, modules, 
    discovery, test, FALSE
  )

  # Now we will order the genes ourselves to prevent duplicate calls to 
  # networkProperties, which can be quite slow.
  if (orderGenesBy == "discovery" && discovery != test)  {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[discovery]]), FALSE, TRUE)
    # Ordering genes by the discovery network however means we have to calculate
    # The network properties in the discovery network
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery, FALSE, orderModules
    ) 
  } else if (orderGenesBy == "none") {
    moduleOrder <- seq_along(props)
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(props[[mi]]$connectivity)
    }
  } else {
    if (missing(orderModules))
      orderModules <- TRUE
    # order modules
    moduleOrder <- 1
    if (length(props) > 1 && orderModules) {
      # Create a matrix of summary expression profiles to measure the similarity
      seps <- matrix(
        0, ncol=length(props), nrow=length(props[[1]]$summaryExpression)
      )
      colnames(seps) <- names(props)
      for (mi in seq_along(props)) {
        seps[,mi] <- props[[mi]]$summaryExpression
      }
      moduleOrder <- hclust(as.dist(1-cor(seps)))$order
    } else {
      moduleOrder <- seq_along(props)
    }
    
    # order genes
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(sort(
        props[[mi]]$connectivity, decreasing=TRUE, na.last=TRUE
      ))
    }
  }
  
  #-----------------------------------------------------------------------------
  # Plot the Module Membership
  #-----------------------------------------------------------------------------
  # now build the Module Membership vector
  MM <- foreach(mi = seq_along(props), .combine=c) %do% {
    props[[mi]]$moduleMembership
  }
  MM <- MM[geneOrder]
  
  # Plot bar chart
  plotBar(
    MM, c(-1,1), moduleAssignments[[discovery]][geneOrder],
    ifelse(MM > 0, palette[2], palette[1]), drawBorders=drawBorders,
    xaxt=plotGeneNames, plotModuleNames=plotModuleNames, 
    xaxt.line=gaxt.line, maxt.line=maxt.line, main=main,
    ylab="Module membership"
  )
}

#' \code{plotConnectivity:} Plot a bar chart of the normalised intramodular 
#' connectivity (see details) for one or more network modules in their discovery
#' dataset, or an independent test dataset.
#' 
#' @rdname plotTopology
#' @export
plotConnectivity <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, orderGenesBy="discovery", orderModules=TRUE,
  plotGeneNames=TRUE, plotModuleNames, main="", 
  palette="#feb24c", drawBorders=FALSE, gaxt.line=-0.5, maxt.line=3, 
  cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderGenesBy <- orderByArgs[pmatch(orderGenesBy, orderByArgs, nomatch=3)]
  
  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- paste0(".temp-objects", getUUID())
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  geneExpression <- unifyDS(dynamicMatLoad(geneExpression, backingpath=tmp.dir))
  coexpression <- unifyDS(dynamicMatLoad(coexpression, backingpath=tmp.dir))
  adjacency <- unifyDS(dynamicMatLoad(adjacency, backingpath=tmp.dir))
  
  # Format optional input data so it doesn't cause cascading error crashes
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  #-----------------------------------------------------------------------------
  # Get ordering of genes in the 'test' dataset by the dataset specified in 
  # 'orderGenesBy'.
  #-----------------------------------------------------------------------------
  # Get the module membership for each module in the test network.
  props <- networkProperties(
    geneExpression, coexpression, adjacency, moduleAssignments, modules, 
    discovery, test, FALSE
  )
  
  # Now we will order the genes ourselves to prevent duplicate calls to 
  # networkProperties, which can be quite slow.
  if (orderGenesBy == "discovery" && discovery != test)  {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[discovery]]), FALSE, TRUE)
    # Ordering genes by the discovery network however means we have to calculate
    # The network properties in the discovery network
    geneOrder <- geneOrder(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery, FALSE, orderModules
    ) 
  } else if (orderGenesBy == "none") {
    moduleOrder <- seq_along(props)
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(props[[mi]]$connectivity)
    }
  } else {
    if (missing(orderModules))
      orderModules <- TRUE
    # order modules
    moduleOrder <- 1
    if (length(props) > 1 && orderModules && !is.null(geneExpression[[test]])) {
      # Create a matrix of summary expression profiles to measure the similarity
      seps <- matrix(
        0, ncol=length(props), nrow=length(props[[1]]$summaryExpression)
      )
      colnames(seps) <- names(props)
      for (mi in seq_along(props)) {
        seps[,mi] <- props[[mi]]$summaryExpression
      }
      moduleOrder <- hclust(as.dist(1-cor(seps)))$order
    } else {
      moduleOrder <- seq_along(props)
    }
    
    # order genes
    geneOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
      names(sort(
        props[[mi]]$connectivity, decreasing=TRUE, na.last=TRUE
      ))
    }
  }
  
  #-----------------------------------------------------------------------------
  # Plot the Connectivity
  #-----------------------------------------------------------------------------
  # now build the (Normalised) Intramodular Connectivity vector
  kIM <- foreach(mi = seq_along(props), .combine=c) %do% {
    # Normalise the connectivity by the maximum. The value has no meaning,
    # just the relative sizes and ranks
    props[[mi]]$connectivity/max(na.omit(props[[mi]]$connectivity))
  }
  kIM <- kIM[geneOrder]
  
  # Plot bar chart
  plotBar(
    kIM, c(0,1), moduleAssignments[[discovery]][geneOrder],
    palette, drawBorders=drawBorders,
    xaxt=plotGeneNames, plotModuleNames=plotModuleNames, 
    xaxt.line=gaxt.line, maxt.line=maxt.line, main=main,
    ylab="Normalised connectivity"
  )
}

#' \code{plotSummaryExpression:} Plot bar charts of the summary expression
#' profiles of one or more network modules in their discovery dataset, or an
#' independent test dataset.
#' 
#' @rdname plotTopology
#' @export
plotSummaryExpression <- function(
  geneExpression, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, orderSamplesBy="test", orderGenesBy="discovery",
  orderModules, plotSampleNames=TRUE, plotModuleNames, 
  main="", palette=c("#762a83", "#1b7837"), drawBorders=FALSE, 
  saxt.line=-0.5, maxt.line=0, cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
  if (is.null(geneExpression))
    stop("Cannot plot gene expression without gene expression data")
  
  if (class(main) != "character")
    stop("'main' must be a characer vector")
  
  orderByArgs <- c("discovery", "test", "none")
  orderGenesBy <- orderByArgs[pmatch(orderGenesBy, orderByArgs, nomatch=3)]
  orderSamplesBy <- orderByArgs[pmatch(orderSamplesBy, orderByArgs, nomatch=3)]
  
  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- paste0(".temp-objects", getUUID())
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  geneExpression <- unifyDS(dynamicMatLoad(geneExpression, backingpath=tmp.dir))
  coexpression <- unifyDS(dynamicMatLoad(coexpression, backingpath=tmp.dir))
  adjacency <- unifyDS(dynamicMatLoad(adjacency, backingpath=tmp.dir))
  
  # Format optional input data so it doesn't cause cascading error crashes 
  geneExpression <- formatGeneExpression(
    geneExpression, length(coexpression), names(coexpression)
  )
  
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  # Sanity check input for consistency.
  checkSets(
    geneExpression, coexpression, adjacency, moduleAssignments, discovery, test
  )
  
  if (is.null(geneExpression[[test]]))
    stop("Cannot plot gene expression without gene expression data")
  
  if (missing(plotModuleNames))
    plotModuleNames <- !missing(modules) && length(modules) > 1
  
  #-----------------------------------------------------------------------------
  # Get ordering of genes and samples in the 'test' dataset by the dataset 
  # specified in 'orderGenesBy' and 'orderSamplesBy'.
  #-----------------------------------------------------------------------------
  props <- networkProperties(
    geneExpression, coexpression, adjacency, moduleAssignments, modules, 
    discovery, test, FALSE
  )
  
  # Determine gene ordering, then sample ordering.
  if (orderGenesBy == "discovery" && discovery != test) {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[discovery]]), FALSE, TRUE)
    
    propsDisc <- networkProperties(
      geneExpression, coexpression, adjacency, moduleAssignments, modules,
      discovery, discovery, FALSE
    ) 
    
    moduleOrder <- names(propsDisc)
    if (length(propsDisc) > 1 && orderModules) {
      # Create a matrix of summary expression profiles to measure the similarity
      seps <- matrix(
        0, ncol=length(propsDisc), 
        nrow=length(propsDisc[[1]]$summaryExpression)
      )
      colnames(seps) <- names(propsDisc)
      for (mi in seq_along(propsDisc)) {
        seps[,mi] <- propsDisc[[mi]]$summaryExpression
      }
      moduleOrder <- names(propsDisc)[hclust(as.dist(1-cor(seps)))$order]
    }
  } else if (orderGenesBy == "none") {
    moduleOrder <- names(props)
  } else {
    if (missing(orderModules))
      orderModules <- ifelse(is.null(geneExpression[[test]]), FALSE, TRUE)
    
    # Order modules and samples by the test network
    moduleOrder <- names(props)
    if (length(props) > 1 && orderModules) {
      seps <- matrix(
        0, ncol=length(props), 
        nrow=length(props[[1]]$summaryExpression)
      )
      colnames(seps) <- names(props)
      for (mi in seq_along(props)) {
        seps[,mi] <- props[[mi]]$summaryExpression
      }
      moduleOrder <- names(props)[hclust(as.dist(1-cor(seps)))$order]
    }
  }
  
  # If test == discovery, we order samples and modules by discovery, and plot
  # the discovery. So we need to make sure the gene expression exists in the
  # discovery
  if (orderSamplesBy == "discovery" && discovery != test) {
    if (is.null(geneExpression[[discovery]])) {
      stop(
        "Expecting gene expression data for the discovery dataset in order",
        " to sort samples"
      )
    }
    if (!exists("propsDisc")) {
      propsDisc <- networkProperties(
        geneExpression, coexpression, adjacency, moduleAssignments, modules,
        discovery, discovery, FALSE
      )
    }
    sampleOrder <- names(sort(
      propsDisc[[1]]$summaryExpression, decreasing=TRUE
    ))
  } else if (orderSamplesBy == "none") {
    sampleOrder <- rownames(geneExpression[[test]])
  } else {
    sampleOrder <- names(sort(
      props[moduleOrder][[1]]$summaryExpression, decreasing=TRUE
    ))
  }
  
  #-----------------------------------------------------------------------------
  # Identify genes and samples from the 'discovery' dataset not present in the 
  # 'test' dataset.
  #-----------------------------------------------------------------------------
  if (all(sampleOrder %nin% rownames(geneExpression[[test]]))) {
    stop(
      "No samples from the 'orderSamplesBy' dataset are present in the",
      " 'test' dataset"
    )
  }
  na.pos <- which(sampleOrder %nin% rownames(geneExpression[[test]]))
  if (length(na.pos) > 0) {
    presentSamples <- sampleOrder[-na.pos]
  } else {
    presentSamples <- sampleOrder
  }
  
  
  #-----------------------------------------------------------------------------
  # Plot the summary expression profiles 
  #-----------------------------------------------------------------------------
  SEP <- foreach(mi = moduleOrder, .combine=cbind) %do% {
    matrix(
      insert.nas(props[[mi]]$summaryExpression[presentSamples], na.pos),
      ncol=1
    )
  }
  colnames(SEP) <- moduleOrder
  rownames(SEP) <- sampleOrder
  
  # Now build the colors
  cols <- matrix(palette[1], nrow(SEP), ncol(SEP))
  cols[SEP > 0] <- palette[2]
  
  # Plot bar chart
  plotMultiBar(
    SEP, rep(list(range(SEP, na.rm=TRUE)), ncol(SEP)),
    cols=cols, drawBorders=drawBorders, main=main, yaxt=plotSampleNames,
    plotModuleNames=plotModuleNames, yaxt.line=saxt.line, maxt.line=maxt.line,
    xlab="Summary Expression"
  )
}

#' \code{plotExpressionLegend:} Plot a legend for the gene expression heatmap
#' for one or more network modules in their discovery dataset, or an independent
#' test dataset.
#' 
#' @rdname plotTopology
#' @export
plotExpressionLegend <- function(
  geneExpression, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, palette=expression.palette(), horizontal=TRUE, 
  legend.main="Expression", legend.tick.size=0.03, laxt.line=2.5, 
  cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input
  #-----------------------------------------------------------------------------
  if (is.null(geneExpression))
    stop("Cannot plot expression legend without gene expression data")
  
  # Temporary directory to store new bigMatrix objects in
  tmp.dir <- paste0(".temp-objects", getUUID())
  dir.create(tmp.dir, showWarnings=FALSE)
  on.exit({
    unlink(tmp.dir, recursive=TRUE)
  }, add=TRUE)
  
  # Unify data structures and load in matrices
  geneExpression <- unifyDS(dynamicMatLoad(geneExpression, backingpath=tmp.dir))
  coexpression <- unifyDS(dynamicMatLoad(coexpression, backingpath=tmp.dir))
  adjacency <- unifyDS(dynamicMatLoad(adjacency, backingpath=tmp.dir))
  
  # If module discovery has not been performed for all datasets, it may be
  # easier for the user to provide a simplified list structuren
  if (!missing(moduleAssignments) && missing(modules)) {
    modules <- unique(moduleAssignments[[discovery]])
  } else if (missing(moduleAssignments) && missing(modules)) {
    modules <- "1"
  } else if (missing(moduleAssignments) && !missing(modules)) {
    stop("'modules' provided but not 'moduleAssignments'")
  }
  
  # Format optional input data so it doesn't cause cascading error crashes
  moduleAssignments <- formatModuleAssignments(
    moduleAssignments, discovery, length(coexpression), names(coexpression),
    ncol(coexpression[[discovery]]), colnames(coexpression[[discovery]])
  )
  
  if (is.null(geneExpression[[test]]))
    stop("Cannot plot summary expression without gene expression data")
  
  #-----------------------------------------------------------------------------
  # Get the range of gene expression for the modules in the test dataset 
  #-----------------------------------------------------------------------------
  modGenes <- getGenes(moduleAssignments, modules, discovery)
  modGenes <- modGenes %sub_in% colnames(geneExpression[[test]])
  if (length(modGenes) == 0)
    stop("None of the module genes are present in the test dataset")
  rg <- range(geneExpression[[test]][,modGenes])
  
  #-----------------------------------------------------------------------------
  # Plot the legend
  #-----------------------------------------------------------------------------
  emptyPlot(c(0,1), c(0,1), bty="n")
  if (all(rg < 0)) {
    addGradientLegend(
      head(palette, length(palette)/2), rg, rg, horizontal, legend.main, 
      xlim=c(0,1), ylim=c(0,1), tick.size=legend.tick.size, axis.line=laxt.line
    )
  } else if (all(rg > 0)) {
    addGradientLegend(
      tail(palette, length(palette)/2), rg, rg, horizontal, legend.main, 
      xlim=c(0,1), ylim=c(0,1), tick.size=legend.tick.size, axis.line=laxt.line
    )
  } else {
    plim <- c(-max(abs(rg)), max(abs(rg)))
    addGradientLegend(
      palette, plim, rg, horizontal, legend.main, xlim=c(0,1), ylim=c(0,1), 
      tick.size=legend.tick.size, axis.line=laxt.line
    )
  }
}

#' \code{plotCoexpressionLegend:} Plot a legend for the gene coexpression.
#' 
#' @rdname plotTopology
#' @export
plotCoexpressionLegend <- function(
  palette=coexpression.palette(),  horizontal=TRUE, legend.main="Coexpression",
  legend.tick.size=0.03, laxt.line=2.5, cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Render legend
  #-----------------------------------------------------------------------------
  emptyPlot(c(0,1), c(0,1), bty="n")
  addGradientLegend(
    palette, c(-1,1), c(-1,1), horizontal, legend.main,
    xlim=c(0,1), ylim=c(0,1), tick.size=legend.tick.size,
    axis.line=laxt.line
  )
}

#' \code{plotAdjacencyLegend:} Plot a legend for the gene adjacencies.
#' 
#' @rdname plotTopology
#' @export
plotAdjacencyLegend <- function(
  palette=adjacency.palette(),  horizontal=TRUE, legend.main="Adjacency",
  legend.tick.size=0.03, laxt.line=2.5, cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
  })
  
  #-----------------------------------------------------------------------------
  # Render legend
  #-----------------------------------------------------------------------------
  emptyPlot(c(0,1), c(0,1), bty="n")
  addGradientLegend(
    palette, c(0,1), c(0,1), horizontal, legend.main,
    xlim=c(0,1), ylim=c(0,1), tick.size=legend.tick.size,
    axis.line=laxt.line
  )
}
