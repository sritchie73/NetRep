
#' Order genes and modules in a network
#' 
#' Clusters network modules by summary expression profile similarity, and 
#' order genes within each subset by connectivity.
#' 
#' @param adjacency Matrix of node adjacencies
#' @param module.labels Optional vector of module assignments for each gene.
#' @param summary.exp Optional \code{data.frame} of summary expression profiles
#'  for each network subset. 
#'  
#' @return a vector of ordered nodes. 
#'
#' @export
orderNetwork <- function(adjacency, module.labels=NULL, summary.exp=NULL) {
  adjacency <- dynamicMatLoad(adjacency)
  if (nrow(adjacency) != ncol(adjacency))
    stop("expecting a square adjacency!")
  
  nGenes <- nrow(adjacency)
  order <- NULL
  
  
  if (is.null(module.labels)) {
    module.labels = rep("a", nGenes)
    module.order = "a"
  } else {
    if (!is.null(summary.exp)) {
      h <- hclust(
        as.dist(
          1 - abs(cor(summary.exp, use="pairwise.complete.obs"))
        )
      )
      module.order <- h$label[h$order]
    } else {
      module.order <- unique(module.labels)
      tryCatch({
        # attempt to order numerically if possible
        module.order <- module.order[order(as.integer(module.order))]
      }, warning = function(w) {
        # leave order as is
      })
    }
  }
  
  for (module in module.order) {
    module.nodes <- which(module.labels == module)
    
    # Order genes within each module by their connectivity to all other genes
    kIM <- adjProps(adjacency, module.nodes)$kIM
    order <- c(order, module.nodes[order(kIM, decreasing=TRUE)])
  }   
  
  return(order)
}

#' Plot Subset legend
#' 
#' Render a horizontal bar that demarcates more clearly where nodes fall into
#' separate modules.
#' 
#' @param module.labels module labels for the nodes in the corresponding 
#'  coexpression heatmap.
#'
#' @export
plotModuleLegend <- function(module.labels, module.colors) {
  if (missing(module.colors)) {
    module.colors = seq_along(module.labels)
  }
    
  nGenes <- length(module.labels)
  
  # create empty plot window
  xlim <- getLimsForManualUsr(c(0, nGenes))
  ylim <- getLimsForManualUsr(c(0, 1))
  plot(
    0, type='n', xaxt='n', yaxt='n', xlim=xlim, ylim=ylim, xlab="", ylab="", 
    frame.plot = FALSE
  )
  
  breaks <- getModuleBreaks(module.labels)
  mods <- rle(module.labels)$values
  
  for (mm in 2:length(breaks)) {
    col = module.colors[match(mods[mm - 1], unique(module.labels))]
    rect(breaks[mm - 1], 0, breaks[mm], 1, col=col, border="black")
  }
  box()
  
}

#' Render gene names for coexpression plot
#' 
#' Create panel containing the axis with just the gene names. Adjust bot.mar to
#' get output correct.
#' 
#' @param gene.labels labels for each gene to render
#' @param bot.mar bottom margin
#' @param cex.axis cex for the axis text
#' 
#' @export
plotGeneNames <- function(gene.labels, bot.mar=5, cex.axis=1) {
  old.mar<- par("mar")
  par(mar=c(bot.mar, old.mar[2], 0, old.mar[4])) 
  
  nGenes <- length(gene.labels)
  xlim <- getLimsForManualUsr(c(0, nGenes))
  ylim <- getLimsForManualUsr(c(0, 0))
  plot(
    0, type='n', xaxt='n', yaxt='n', xlim=xlim, ylim=ylim, xlab="", ylab="", 
    frame.plot = FALSE
  )
  axis(at=(1:nGenes)-0.5, side=1, labels=gene.labels, las=2, cex.axis=cex.axis)
  par(mar=old.mar)
}

#' Plot legend for coexpression
#' 
#' @param heatmap.gradient A vector of colors that was used in 
#'  \code{\link{plotCoexpression}}
#' @param nBins number of bins to show on the legend.
#' @param cex.title cex for the title text
#' @param cex.legend cex for the axis text
#' 
#' @export
plotCoexpLegend <- function(
  heatmap.gradient, nBins=11, title="Pearson\nCorrelation", cex.title=1.4, 
  cex.legend=1
) {
  if(missing(heatmap.gradient)) {
    heatmap.gradient <- custom.palette()
  }
  grad <- colorRampPalette(heatmap.gradient)(nBins)
  
  old.mar <- par('mar')
  
  par(mar=c(1,4,4,1))
  
  xlim <- getLimsForManualUsr(c(0, 1))
  ylim <- getLimsForManualUsr(c(0, nBins))
  plot(
    0, type='n', xaxt='n', yaxt='n', xlim=xlim, ylim=ylim, xlab="", ylab="", 
    frame.plot = FALSE
  )
  
  ltext <- format(seq(-1, 1, length=nBins+1), digits=2)
  
  for (ii in seq_along(grad)) {
    rect(
      0, ii-1, 1, ii, border=NA, col=grad[ii]  
    )
  }
  box()
  axis(side=2, at=c(0:nBins), labels=ltext, las=2, cex.axis=cex.legend)
  mtext(text=title, side=3, cex=cex.title, line=1)
  
  par(mar=old.mar)
}

#' Coexpression plot
#' 
#' Plot a triangle of the coexpression
#' 
#' @param coexpresion a square matrix containing the pairwise gene coexpression.
#' @param module.labels an optional named vector assigning each node to a module.
#' @param heatmap.gradient A vector of colors to use for (or interpolate over)
#'        to render the coexpression.
#' @export
plotCoexpression <- function(
  coexpression, module.labels=NULL, heatmap.gradient
) {
  coexpression <- dynamicMatLoad(coexpression)
  if (nrow(coexpression) != ncol(coexpression))
    stop("expecting a square edge.matrix!")
  
  nGenes <- nrow(coexpression)
  
  if(missing(heatmap.gradient)) {
    heatmap.gradient <- custom.palette()
  }
  nColBins <- 255 
  colGrad <- colorRampPalette(heatmap.gradient)(nColBins)
  edgeBins <- seq(-1, 1, length=nColBins+1)
  
  # create empty plot window
  xlim <- getLimsForManualUsr(c(0, nGenes))
  ylim <- getLimsForManualUsr(c(0, nGenes/2))
  plot(
    0, type='n', xaxt='n', yaxt='n', xlim=xlim, ylim=ylim, xlab="", ylab="", 
    frame.plot = FALSE
  )
  
  # render correlation squares / triangles
  for (ii in 1:nGenes) {
    for (jj in 1:ii) {
      plotRow <- ii - jj
      topy <- (plotRow + 1)/2
      # If we're on the diagonal, plot a triangle, otherwise a diamond
      if (plotRow == 0) {
        boty <- 0
      } else {
        boty <- topy - 1
      }
      xOffset <- plotRow/2
      rightx <- jj + xOffset
      leftx <- rightx - 1
      
      polygon(
        x=c(leftx, leftx+0.5, rightx, leftx+0.5, leftx),
        y=c(topy-0.5, topy, topy-0.5, boty, topy-0.5),
        col=findColInGrad(coexpression[ii, jj], edgeBins, colGrad),
        border=NA
      )
    }
  }
  
  # render module boundaries
  if (!is.null(module.labels)) {
    breaks <- getModuleBreaks(module.labels)
    for (mm in 2:length(breaks)) {
      mheight <- (breaks[mm] - breaks[mm - 1])/2
      halfway <-  mheight + breaks[mm - 1]
      polygon(
        x=c(breaks[mm - 1], breaks[mm], halfway, breaks[mm - 1]),
        y=c(0, 0, mheight, 0)
      )
    }
  }
  
  # render border of plot
  polygon(x=c(0, nGenes, nGenes/2, 0), y=c(0, 0, nGenes/2, 0))
  
}

# Get edge weight color from gradient.
# 
# Expects `edgeBins` to be in ascending order, and `edgeBins` to be 1 longer 
# than `colors`.
#
findColInGrad <- function(weight, edgeBins, colors) {  
  # Handle cases where edge weight is outside range
  if (weight <= edgeBins[1]){
    return(colors[1])
  }
  if (weight >= tail(edgeBins, 1)) {
    return(tail(colors, 1))
  }
  # Search for bin weight falls into
  for (i in 1:(length(edgeBins) - 1)) {
    if (weight >= edgeBins[i] & weight <= edgeBins[i + 1]) {
      return(colors[i])
    }
  }
}

#' Modified Brewer Palette
#' 
#' RColorBrewer palette "RdYlBu" with the middle color replaced with white. This
#' gives a nicer contrast than the "RdBu" palette
#' 
custom.palette <- function() {
  c(
    "#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFFF", "#FEE090", 
    "#FDAE61", "#F46D43", "#D73027", "#A50026"
   )
}

#' Get the plot limits to set for the desired plot window
#'
#' \code{plot} manually adds an extra region to the plot on top of the given
#' 'xlim' and 'ylim', amounting to 4% of the plot region in either direction.
#' This function tells you what limits to set so that the boundaries of the plot
#' are precisely the min and max you want.
#' 
#' @param dlim the limits you want to set
#' @return
#'  the min and max to set for the desired plot limits  
getLimsForManualUsr <- function(dlim) {
  A = matrix(c(1.04, -0.04, -0.04, 1.04), nrow=2, ncol=2)
  B = matrix(dlim, nrow=2, ncol=1)
  as.vector(solve(A, B))
}

# Identifies the break points between modules on the x axis.
getModuleBreaks <- function(module.labels) {
  breaks <- rle(module.labels)
  break.points <- sapply(1:length(breaks$length), function(i) {
    sum(breaks$length[1:i])
  })  
  break.points <- c(0, break.points)
  break.points
}
