#' Plot the topology of a gene coexpression module
#' 
#' @description 
#' Plot the coexpression, adjacency, (normalised) intramodular connectivity, 
#' module membership, gene expression, and summary expression profiles of one or
#' more network modules in the discovery dataset or independent test dataset.
#' 
#' Individual components of the module plot can be plotted using 
#' \code{\link{plotCoexpression}}, \code{\link{plotAdjacency}}, 
#' \code{\link{plotConnectivity}}, \code{\link{plotModuleMembership}}, 
#' \code{\link{plotExpression}}, and \code{\link{plotSummaryExpression}}.
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
#' @param orderGenesBy one of "discovery", "test", or "none" indicating which
#'  dataset the function should use to order the genes (see details).
#' @param orderSamplesBy one of "discovery", "test", or "none" indicating which
#'  dataset the function should use to order the samples (see details).
#' @param orderModules logical; if \code{TRUE} modules will be ordered by the 
#'  similarity of their summary expression profiles (see details). The default
#'  is \code{TRUE} if gene expression for the \code{test} dataset is provided.
#' @param plotGeneNames logical; if \code{TRUE}, plot the gene names on the 
#'  bottom axis.
#' @param plotSampleNames logical; if \code{TRUE}, plot the sample names on the
#'  left axis next to the summary expression profiles.
#' @param plotModuleNames logical; if \code{TRUE}, plot the module names on the
#'  bottom axis and above the summary expression profiles. By default, the 
#'  module names are rendered only if multiple \code{modules} are specified.
#' @param main title for the plot.
#' @param drawBorders logical; if \code{TRUE}, borders are drawn around the 
#'  connectivity, module membership, and summary expression bar plots.
#' @param gaxt.line the number of lines into the bottom margin at which the gene
#'  names will be drawn.
#' @param saxt.line the number of lines into the left margin at which the sample
#'  names will be drawn.
#' @param maxt.line the number of lines into the bottom margin at which the 
#'  module names will be drawn.
#' @param legend.tick.size size of the ticks on each axis legend relative to the
#'  size of the coexpression, adjacency, and gene expression heatmaps.
#' @param laxt.line the distance from the legend to render the legend axis 
#'  labels, as multiple of \code{legend.tick.size}.
#' @param cex.axis relative size of the gene and sample names.
#' @param cex.lab relative size of the module names, legend titles, and axis
#'  labels.
#' @param cex.main relative size of the plot title.
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
#'   By default, samples are ordered in descending order of the summary 
#'   expression profile for the left-most module appearing on the plot (see 
#'   \code{\link{sampleOrder}}. By default, the summary expression profile is
#'   calculated in the \code{test} dataset. This behaviour can be changed
#'   through the \code{orderSamplesBy} argument, however setting
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
#' \subsection{Plot customisation:}{
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
#'   \code{cex.main} controls the relative text size of the plot title
#'   (specified by the \code{main} argument). \code{cex.axis} controls the
#'   relative text size of the gene and sample labels. \code{cex.lab} controls
#'   the relative text size of the bar plot axis labels, module labels, and the
#'   legend titles.
#'   
#'   The rendering of gene, sample, and module names can be disabled by setting
#'   \code{plotGeneNames}, \code{plotSampleNames}, and \code{plotModuleNames} to
#'   \code{FALSE}.
#'   
#'   The \code{drawBorders} argument controls whether borders are drawn around
#'   the connectivity, module membership, or summary expression bar plots.
#' }
#' 
#' @seealso
#' \code{\link{plotCoexpression}} 
#' \code{\link{plotAdjacency}}
#' \code{\link{plotConnectivity}}
#' \code{\link{plotModuleMembership}}
#' \code{\link{plotExpression}}
#' \code{\link{plotSummaryExpression}}
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
#' plotModule(geA, coexpA, adjA, moduleAssignments, modules="2")
#' 
#' ## Example 2: Plot an arbitrary set of genes in cohort A
#' plotModule(geA[,1:10], coexpA[1:10, 1:10], adjA[1:10, 1:10])
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
#' plotModule(
#'   geneExpression, coexpression, adjacency, moduleAssignments,
#'   modules=c("3", "7"), discovery="adipose", test="liver"
#' )
#' 
#' # clean up bigMatrix files from examples
#' unlink("*_bm*")
#' }
#' 
#' @name plotModule
#' @export
plotModule <- function(
  geneExpression, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1, orderSamplesBy="test", orderGenesBy="discovery",
  orderModules, plotGeneNames=TRUE, plotSampleNames=TRUE, plotModuleNames,
  main="Module Topology", drawBorders=FALSE, gaxt.line=-0.5, 
  saxt.line=-0.5, maxt.line, legend.tick.size=0.04, 
  laxt.line=2.5, cex.axis=0.8, cex.lab=1, cex.main=1.2
) {
  #-----------------------------------------------------------------------------
  # Set graphical parameters
  #-----------------------------------------------------------------------------
  old.par <- par(c("cex.axis", "cex.lab", "cex.main", "mar", "oma"))
  par(cex.axis=cex.axis)
  par(cex.lab=cex.lab)
  par(cex.main=cex.main)
  # make sure to restore old values once finishing the plot
  on.exit({
    par(cex.axis=old.par[[1]])
    par(cex.lab=old.par[[2]])
    par(cex.main=old.par[[3]])
    par(mar=old.par[[4]])
    par(oma=old.par[[5]])
    layout(1)
  })
  
  #-----------------------------------------------------------------------------
  # Validate user input and unify data structures
  #-----------------------------------------------------------------------------
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
      if (!is.null(geneExpression[[discovery]])) {
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
      } else {
        warning(
          "Gene expression missing from the 'discovery' dataset, modules will ",
          "be ordered as provided"
        )
      }
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
      if (!is.null(geneExpression[[test]])) {
        seps <- matrix(
          0, ncol=length(props), 
          nrow=length(props[[1]]$summaryExpression)
        )
        colnames(seps) <- names(props)
        for (mi in seq_along(props)) {
          seps[,mi] <- props[[mi]]$summaryExpression
        }
        moduleOrder <- names(props)[hclust(as.dist(1-cor(seps)))$order]
      } else {
        warning(
          "Gene expression missing from the 'test' dataset, modules will be ",
          "ordered as provided"
        )
      }
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
  } else if (!is.null(geneExpression[[test]])) {
    if (orderSamplesBy == "none") {
      sampleOrder <- rownames(geneExpression[[test]])
    } else {
      sampleOrder <- names(sort(
        props[moduleOrder][[1]]$summaryExpression, decreasing=TRUE
      ))
    }
  }
  
  #-----------------------------------------------------------------------------
  # Identify genes and samples from the 'discovery' dataset not present in the 
  # 'test' dataset.
  #-----------------------------------------------------------------------------
  if (!is.null(geneExpression[[test]])) {
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
  }

  # Handle genes not present in the test dataset
  na.pos.x <- which(geneOrder %nin% colnames(coexpression[[test]]))
  if (length(na.pos.x) > 0) {
    presentGenes <- geneOrder[-na.pos.x]
  } else {
    presentGenes <- geneOrder
  }
  
  #-----------------------------------------------------------------------------
  # Set up other property vectors and datasets
  #-----------------------------------------------------------------------------
  # (Normalised) Intramodular Connectivity vector
  kIM <- foreach(mi = seq_along(props), .combine=c) %do% {
    # Normalise the connectivity by the maximum. The value has no meaning,
    # just the relative sizes and ranks
    props[[mi]]$connectivity/max(na.omit(props[[mi]]$connectivity))
  }
  kIM <- kIM[geneOrder]
  
  if (!is.null(geneExpression[[test]])) {
    # now build the Module Membership vector
    MM <- foreach(mi = seq_along(props), .combine=c) %do% {
      props[[mi]]$moduleMembership
    }
    MM <- MM[geneOrder]
    MM.cols <- rep(coexpression.palette()[1], length(MM))
    MM.cols[MM > 0] <- tail(coexpression.palette(), 1)
    
    # Summary Expression profiles matrix
    SEP <- foreach(mi = moduleOrder, .combine=cbind) %do% {
      matrix(
        insert.nas(props[[mi]]$summaryExpression[presentSamples], na.pos.y),
        ncol=1
      )
    }
    colnames(SEP) <- moduleOrder
    rownames(SEP) <- sampleOrder
    
    # Now build the colors
    cols <- matrix(expression.palette()[1], nrow(SEP), ncol(SEP))
    cols[SEP > 0] <- tail(expression.palette(), 1)
  }
  
  #-----------------------------------------------------------------------------
  # Set up plotting region
  #-----------------------------------------------------------------------------
  gaxt <- NULL
  if (plotGeneNames)
    gaxt <- geneOrder
  
  if (missing(maxt.line) && !plotGeneNames) {
    maxt.line <- -0.5
  } else if (missing(maxt.line) && plotGeneNames) {
    maxt.line <- 3
  }
  
  if (is.null(geneExpression[[test]])) {
    # set up plot layout
    layout(
      mat=matrix(1:3, ncol=1), heights=c(0.4, 0.4, 0.2)
    )
  } else {
    # set up plot layout
    SEPw <- min(0.2 + (length(modules) - 1)* 0.1, 0.5) 
    layout(
      mat=matrix(c(rep(8, 5), 7, 1:6), ncol=2), 
      heights=c(0.7/3, 0.7/3, 0.12, 0.12, 0.06, 0.7/3),
      widths=c(SEPw, 1-SEPw)
    )
  }
  par(oma=old.par[["mar"]]+old.par[["oma"]])
  
  #-----------------------------------------------------------------------------
  # Plot topology components
  #-----------------------------------------------------------------------------
  # Plot coexpression
  par(mar=c(1, 1, 1, 1))
  plotTriangleHeatmap(
    coexpression[[test]][presentGenes , presentGenes], 
    coexpression.palette(), c(-1, 1), 
    moduleAssignments[[discovery]][geneOrder], na.pos.x, plotLegend=TRUE, 
    main="", legend.main="Coexpression", plotModuleNames=FALSE, 
    legend.tick.size=legend.tick.size, laxt.line=laxt.line, 
    legend.line=0.1, maxt.line=maxt.line
  )
  mtext(main, side=3, line=1, cex=par('cex.main'), font=2, xpd=NA)
  
  # Plot adjacency
  par(mar=c(1, 1, 1, 1))
  plotTriangleHeatmap(
    adjacency[[test]][presentGenes , presentGenes], adjacency.palette(), 
    c(0, 1), moduleAssignments[[discovery]][geneOrder], na.pos.x, 
    plotLegend=TRUE, main="", legend.main="Adjacency", 
    plotModuleNames=FALSE, legend.tick.size=legend.tick.size, 
    laxt.line=laxt.line, legend.line=0.1, maxt.line=maxt.line
  )
  
  # Plot Intamodular Connectivity
  par(mar=c(1,1,1,1))
  if (is.null(geneExpression[[test]])) {
    plotBar(
      kIM, c(0,1), moduleAssignments[[discovery]][geneOrder], "#feb24c", 
      drawBorders=drawBorders, plotModuleNames=plotModuleNames, 
      xaxt=plotGeneNames, xaxt.line=gaxt.line, main="",
      ylab="Normalised connectivity", maxt.line=maxt.line
    )
  } else {
    plotBar(
      kIM, c(0,1), moduleAssignments[[discovery]][geneOrder], "#feb24c", 
      drawBorders=drawBorders, plotModuleNames=FALSE, main="", xaxt=FALSE,
      ylab="Normalised\nconnectivity", maxt.line=maxt.line
    )
  }
  
  if (!is.null(geneExpression[[test]])) {
    # Plot Module Membership
    par(mar=c(1, 1, 1, 1))
    plotBar(
      MM, c(-1,1), moduleAssignments[[discovery]][geneOrder], MM.cols,
      drawBorders=drawBorders, plotModuleNames=FALSE, main="", xaxt=FALSE,
      ylab="Module\nMembership", maxt.line=maxt.line
    )
    
    # Plot the gene expression
    ge <- geneExpression[[test]][presentSamples, presentGenes]
    range.ge <- range(ge)
    if (all(range.ge > 0)) {
      palette <- tail(expression.palette(), length(expression.palette())/2)
    } else if (all(range.ge < 0)) {
      palette <- head(expression.palette(), length(expression.palette())/2)
    } else {
      palette <- expression.palette()
      range.pal <- c(-max(abs(range.ge)), max(abs(range.ge)))
    }
    
    par(mar=c(1,1,1,1))
    emptyPlot(c(0,1), c(0,1), bty="n")

    yaxt <- NULL
    if (plotSampleNames)
      yaxt <- sampleOrder
    par(mar=c(1, 1, 1, 1))
    plotSquareHeatmap(
      ge, palette, vlim=range.pal, legend.lim=range.ge,
      moduleAssignments[[discovery]][geneOrder], na.pos.x, na.pos.y, 
      xaxt=gaxt, yaxt=NULL, plotLegend=FALSE, main="",
      legend.main="", plotModuleNames=plotModuleNames,
      xaxt.line=gaxt.line, maxt.line=maxt.line
    )
    
    # Plot gene expression legend
    if (all(range.ge < 0)) {
      addGradientLegend(
        head(palette, length(palette)/2), range.ge, range.ge, TRUE, 
        xlim=c(0.5+ncol(ge)*0.1,ncol(ge)+0.5-ncol(ge)*0.1), 
        ylim=c(nrow(ge)+0.5+nrow(ge)*0.2,nrow(ge)+0.5+nrow(ge)*0.3), 
        tick.size=legend.tick.size, 
        main="Gene expression", axis.line=laxt.line
      )
    } else if (all(range.ge > 0)) {
      addGradientLegend(
        tail(palette, length(palette)/2), range.ge, range.ge, TRUE, 
        xlim=c(0.5+ncol(ge)*0.1,ncol(ge)+0.5-ncol(ge)*0.1), 
        ylim=c(nrow(ge)+0.5+nrow(ge)*0.2,nrow(ge)+0.5+nrow(ge)*0.3), 
        tick.size=legend.tick.size, 
        main="Gene expression", axis.line=laxt.line
      )
    } else {
      plim <- c(-max(abs(range.ge)), max(abs(range.ge)))
      addGradientLegend(
        palette, plim, range.ge, TRUE, main="Gene expression",
        xlim=c(0.5+ncol(ge)*0.1,ncol(ge)+0.5-ncol(ge)*0.1), 
        ylim=c(nrow(ge)+0.5+nrow(ge)*0.2,nrow(ge)+0.5+nrow(ge)*0.3),  
        tick.size=legend.tick.size, 
        axis.line=laxt.line
      )
    }
    
    # Plot bar chart
    xlab <- "Summary Expression"
    if (length(modules) == 1) 
      xlab <- gsub(" ", "\n", xlab)
    par(mar=c(1, 1, 1, 1))
    plotMultiBar(
      SEP, rep(list(range(SEP, na.rm=TRUE)), ncol(SEP)),
      cols=cols, drawBorders=drawBorders,
      yaxt=plotSampleNames, plotModuleNames=plotModuleNames, 
      yaxt.line=saxt.line, maxt.line=0, xlab=xlab, 
      cex.modules=par("cex.lab")*0.7
    )
  }
}