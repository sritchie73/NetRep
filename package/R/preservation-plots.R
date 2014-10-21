#' Generate preservation plot components
#' 
#' @details
#' All matrix arguments will take \code{matrix}, \code{big.matrix}, or 
#' \code{big.matrix} descriptors.
#' 
#' @return
#' Outputs files for the following:
#' \itemize{
#' \item coexpression heatmap
#' \item coexpression legend
#' \item intramodular connectivity
#' \item intramodular module membership
#' \item expression heatmap
#' \item summary expression profile
#' }
#' 
#' @param gene.expr gene expression matrix. Expects rows to be samples, and 
#'  columns to be genes. 
#' @param coexpression matrix of pairwise gene coexpression.
#' @param adjacency matrix of pairwise gene adjacencies.
#' @param moduleLabels vector of labels assigning each gene to a module.
#' @param module module to plot preservation for.
#' @param is.relative logical; is the gene expression relative to some other 
#'   measure? If \code{TRUE}, a divergent palette is selected when plotting the
#'   expression, otherwise a sequential palette is chosen.
#' @param cex.axis \code{cex} for the axis tick labels
#' @param cex.title \code{cex} for the titles
#' @param cex.lab \code{cex} for the axis labels
#' @param dpi resolution for each of the panels.#' 
#' @export
preservationPlot <- function(
  gene.expr, coexpression, adjacency, moduleLabels, module,
  is.relative=TRUE, cex.axis=0.8, cex.title=0.9, dpi=100, cex.lab=0.6
) {
  old.mar <- par('mar')
  
  # Load matrices
  ge <- dynamicMatLoad(gene.expr)
  coexp <- dynamicMatLoad(coexpression)
  adj <- dynamicMatLoad(adjacency)
  
  # Get module subset
  modNodes <- which(moduleLabels == module)
  ge <- ge[,modNodes]
  coexp <- coexp[modNodes, modNodes]
  adj <- adj[modNodes, modNodes]
  
  # Order genes in decreasing order of connectivity
  go <- orderNetwork(adj)
  ge <- ge[, go]
  coexp <- coexp[go, go]
  adj <- adj[go, go]

  # Order samples in decreasing order by summary expression profile
  so <- orderSamples(ge)
  ge <- ge[so,]
  
  # cast back to big.matrix
  adj <- dynamicMatLoad(adj)
  coexp <- dynamicMatLoad(coexp)
  ge <- dynamicMatLoad(ge)
  
  nGenes <- length(go)
  
  titleh <- 0.04 # There are 4 of these
  trih <- 0.2 # There are 2 of these
  exph <- 0.2 # There is 1 of these
  barh <- 0.058 # There are 2 of these
  seph <- titleh*1.85 # There is 1 of these
  axh <- titleh*1.25 # There is 1 of these
  
  # debug layout heights sum to 1
  # print(titleh*4 + trih*2 + exph + barh*2 + seph + axh) 
  
  layout(
    matrix(1:33, ncol=3), widths=c(0.2, 0.65, 0.15), 
    heights=c(
      titleh, trih, titleh, trih, titleh, barh, titleh, barh, seph, exph, axh
    )
  )
  on.exit({ dev.off(); par(mar=old.mar) }) # Turn off devices if there's a crash
  
  title.line <- 0
  
  #----------------------------------------------
  # First column of plots, mostly blank or axes
  #----------------------------------------------
  
  par(mar=c(0,0,0,0))
  nullPlot() # empty space next to title above adjacency heatmap
  nullPlot() # empty space next to adjacency heatmap
  nullPlot() # empty space next to title above coexpression heatmap
  nullPlot() # empty space next to coexpression heatmap
  nullPlot() # empty space next to title above intramodular connectivity plot
  nullPlot() # empty space reserved for intramodular connectivity axis
  nullPlot() # empty space next to title above intramodular module membership
  nullPlot() # empty space reserved for intramodular module membership axis
  
  # Plot the summary expression
  par(mar=c(0,0.4,0,1))
  nullPlot() # empty space reserved for title
  plotSummaryExpression(
    ge, is.relative=is.relative, cex.axis=cex.axis, cex.title=cex.title*0.7,
    ylab=""
  )
  
  par(mar=c(0,0,0,0))
  nullPlot() # empty space reserved for summary expression profile plot axis

  
  #----------------------------------------------
  # second column of main plots
  #----------------------------------------------
  
  # Plot adjacency heatmap
  par(mar=c(1,0,0,0))
  nullPlot() # empty space reserved for title
  plotCoexpression(adj, cex.title=cex.title, cex.lab=cex.lab, adjacency=TRUE)
  
  # Plot coexpression heatmap
  nullPlot() # empty space reserved title
  plotCoexpression(coexp, cex.title=cex.title, cex.lab=cex.lab)
  
  # Plot intramodular connectivity
  nullPlot() # empty space reserved for title
  plotKIM(adj, cex.title=cex.title, cex.axis=cex.axis, cex.lab=cex.lab)
  
  # Plot intramodular module membership
  par(mar=c(0,0,0,0))
  nullPlot() # empty space reserved for title
  plotKME(ge, cex.title=cex.title, cex.axis=cex.axis, cex.lab=cex.lab)
  
  # Plot expression matrix
  nullPlot() # empty space reserved for title
  plotExpression(
    ge, is.relative=is.relative, cex.title=cex.title, cex.lab=cex.lab 
  )
  
  # Empty plot underneath the expression (for the SEP axis)
  nullPlot()
  
  #----------------------------------------------
  # third column for legends
  #----------------------------------------------
  
  nullPlot() # empty space next to title above adjacency heatmap
  
  # Plot adjacency legend
  par(mar=c(2,3,2,0.4))
  plotLegend( 
    gradient=custom.palette(FALSE), 
    range=c(0, 1), cex.axis=cex.axis, cex.title=cex.title*0.7, main="", nTicks=4
  )
  
  par(mar=c(0,0,0,0))
  nullPlot() # empty space next to title above coexpression heatmap
  
  # Plot coexpression legend
  par(mar=c(2,3,2,0.4))
  plotLegend( 
    gradient=custom.palette(TRUE), 
    range=c(-1, 1), cex.axis=cex.axis, cex.title=cex.title*0.7, main="", 
    nTicks=5
  )
  
  par(mar=c(0,0,0,0))
  nullPlot() # empty space next to title above intramodular connectivity
  nullPlot() # empty space next to intramodular connectivity
  nullPlot() # empty space next to title above intramodular module membership
  nullPlot() # empty space next to intramodular module membership
  nullPlot() # empty space next to title above expression heatmap
  
  # Plot expression legend
  par(mar=c(2,3,2,0.4))
  plotLegend( 
    gradient=expression.palette(is.relative), 
    range=rangeBigMatrix(ge), cex.axis=cex.axis, cex.title=cex.title*0.7, 
    main="", nTicks=5
  )
  
  par(mar=c(0,0,0,0))
  nullPlot() # empty space next to SEP axis.
  
  par(mar=old.mar)
  on.exit({ }) # Remove dev.off statement now we've succeeded
  
}

#------------------------------------------------------------------
# Functions to organise genes and samples in a visually useful way
#------------------------------------------------------------------

#' Order genes and modules in a network
#' 
#' Clusters network modules by summary expression profile similarity, and 
#' order genes within each subset by connectivity.
#' 
#' @param adjacency Matrix of gene adjacencies
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

#' Order gene expression samples based on the summary expression profile
#' 
#' @param gene.expr matrix of gene expression data. Expects columns to be the 
#'  genes.
#' @return a vector of ordered nodes. 
#'
#' @export
orderSamples <- function(gene.expr) {
  gene.expr <- dynamicMatLoad(gene.expr)
  scaled <- scaleBigMatrix(gene.expr)
  SEP <- dataProps(gene.expr, scaled, 1:ncol(gene.expr))$SEP
  o <- order(SEP, decreasing=TRUE)
  o
}

#------------------------------------------------------------------
# Individual component plots. Also useful by themselves.
#------------------------------------------------------------------

#' Plot a vertical gradient legend
#' 
#' @param gradient A vector of colors to create the gradient
#' @param range range of values the colors fall over
#' @param cex.axis cex for the axis text
#' @param cex.title cex for the title text
#' @param main title for the legend
#' @param nTicks number of tick marks for the axis
#' 
#' @export
plotLegend <- function(
  gradient, range, cex.axis=0.8, cex.title=1, main="Legend", nTicks=9
) {
  if(missing(gradient)) {
    gradient <- custom.palette()
  }
  nColBins <- 255
  grad <- colorRampPalette(gradient)(nColBins)
  
  nullPlot(c(0, 1), range)
  
  axis.locs <- seq(range[1], range[2], length=nTicks)
  axis.text <- format(axis.locs, digits=2)
  
  binLocs <- seq(range[1], range[2], length=nColBins + 1)
  for (ii in 2:length(binLocs)) {
    rect(
      xleft=0, 
      ybottom=binLocs[ii-1], 
      xright=1, 
      ytop=binLocs[ii], 
      border=NA, col=grad[ii-1]  
    )
  }
  box()
  axis(side=2, at=axis.locs, labels=axis.text, las=2, cex.axis=cex.axis)
  mtext(text=main, side=3, cex=cex.title, line=1, adj=1, padj=1)
}

#' Plot the Intramodular Connectivity
#' 
#' @param adjacency Matrix of gene adjacencies
#' @param module.labels Optional vector of module assignments for each gene.
#' @param cex.axis cex for the axis text
#' @param cex.lab cex for the axis label
#' @param cex.title cex for the title text
#' @param main title for the plot
#' 
#' @export
plotKIM <- function(
  adjacency, module.labels, col="orange", cex.axis=1, cex.title=1.4, 
  main="Intramodular Connectivity", xlab="genes", cex.lab=cex.axis
) {
  adjacency <- dynamicMatLoad(adjacency)
  if (nrow(adjacency) != ncol(adjacency))
    stop("expecting a square adjacency!")
  
  nGenes <- nrow(adjacency)
  
  if (missing(module.labels)) {
    module.labels <- rep("a", nGenes)
  }
  
  poke(adjacency)
  kIM <- rep(0, nGenes)
  
  for (mm in unique(module.labels)) {
    modGenes <- which(module.labels == mm)
    kIM[modGenes] <- adjProps(adjacency, modGenes)$kIM
  }
  
  myBarPlot(
    kIM, cols=col, height.lim=c(0, max(kIM)), cex.axis=cex.axis, 
    cex.title=cex.title, main=main, bar.lab=xlab, cex.lab=cex.lab
  )
}

#' Plot the Intramodular Module Membership
#' 
#' @param gene.expr matrix of gene expression data. Expects columns to be the 
#'  genes.
#' @param module.labels Optional vector of module assignments for each gene.
#' @param heatmap.gradient A vector of colors to use for (or interpolate over)
#'        to render the coexpression.
#' @param cex.axis cex for the axis text
#' @param cex.lab cex for the axis label
#' @param cex.title cex for the title text
#' @param main title for the plot
#' @param xlab label for the x axis
#'
#' @export
plotKME <- function(
  gene.expr, module.labels, heatmap.gradient, cex.axis=1, cex.title=1.4,
  main = "Intramodular Module Membership", xlab="genes", cex.lab=cex.axis
) {
  gene.expr <- dynamicMatLoad(gene.expr)
  scaled <- scaleBigMatrix(gene.expr)
  
  nGenes <- ncol(gene.expr)
  
  if (missing(module.labels)) {
    module.labels <- rep("a", nGenes)
  }
  
  if(missing(heatmap.gradient)) {
    heatmap.gradient <- custom.palette()
  }
  nColBins <- 2
  colGrad <- colorRampPalette(heatmap.gradient)(nColBins)
  edgeBins <- seq(-1, 1, length=nColBins+1)
  
  poke(gene.expr)
  poke(scaled)
  kME <- rep(0, nGenes)
  
  for (mm in unique(module.labels)) {
    modGenes <- which(module.labels == mm)
    kME[modGenes] <- dataProps(gene.expr, scaled, modGenes)$kME
  }
  
  cols <- sapply(kME, findColInGrad, edgeBins, colGrad)
  
  myBarPlot(
    kME, cols=cols, cex.axis=cex.axis, height.lim=c(-1, 1), cex.title=cex.title,
    main=main, bar.lab=xlab, cex.lab=cex.lab
  )
  abline(h=0) 
}

#' Plot the Summary Expression Profile of a Module
#' 
#' Vertical barplot.
#' 
#' @param coexpresion a square matrix containing the pairwise gene coexpression.
#' @param module.labels an optional named vector assigning each node to a module.
#' @param heatmap.gradient A vector of colors to use for (or interpolate over)
#'        to render the coexpression.
#' @param is.relative logical; is the expression relative, i.e. centred around 0?
#'   This affects automated color gradient choice.
#' @param cex.axis cex for the axis text
#' @param cex.lab cex for the axis label
#' @param cex.title cex for the title text
#' @param main title for the plot
#' 
#' @export
plotSummaryExpression <- function(
  gene.expr, heatmap.gradient, expression.range, is.relative=TRUE, cex.axis=1,
  cex.title=1.4, main="Summary\nExpression", ylab="samples", cex.lab=cex.axis
) {
  gene.expr <- dynamicMatLoad(gene.expr)
  scaled <- scaleBigMatrix(gene.expr)
  
  if (missing(heatmap.gradient)) {
    heatmap.gradient <- expression.palette(is.relative)
  }
  
  nColBins <- 2
  colGrad <- colorRampPalette(heatmap.gradient)(nColBins)
  exprBins <- seq(-1, 1, length=nColBins+1)
  
  poke(gene.expr)
  poke(scaled)
  SEP <- dataProps(gene.expr, scaled, 1:ncol(gene.expr))$SEP
  
  if (missing(expression.range)) {
    expression.range <- range(SEP)
  }
  
  cols <- sapply(SEP, findColInGrad, exprBins, colGrad)
  
  myBarPlot(
    SEP, cols=cols, cex.axis=cex.axis, height.lim=expression.range, horiz=TRUE, 
    border=NA, bar.width=1, cex.title=cex.title, main=main, bar.lab=ylab, 
    cex.lab=cex.lab
  )
}

#' Plot a heatmap of the gene expression
#' 
#' @param gene.expr matrix of gene expression data. Expects columns to be the 
#'  genes.
#' @param heatmap.gradient A vector of colors to use for (or interpolate over)
#'        to render the gene expression. 
#' @param expression.range range of values the gene expression can take.
#' @param is.relative logical; is the expression relative, i.e. centred around 0?
#'   This affects automated color gradient choice.
#' @param cex.title cex for the title text
#' @param cex.lab cex for the axis label
#' @param main title for the plot   
#' @param xlab label for the x axis
#' @param ylab label for the y axis
#'
#' @importFrom RColorBrewer brewer.pal
#' @export
plotExpression <- function(
  gene.expr, heatmap.gradient, expression.range, is.relative=TRUE, 
  cex.title=1.4, main="Gene Expression", xlab="genes", ylab="samples", 
  cex.lab=1
) {
  gene.expr <- dynamicMatLoad(gene.expr)
  
  if (missing(heatmap.gradient)) {
    heatmap.gradient <- expression.palette(is.relative)
  }
  
  if (missing(expression.range)) {
    poke(gene.expr)
    expression.range <- rangeBigMatrix(gene.expr)
  }
  nColBins <- 255 
  colGrad <- colorRampPalette(heatmap.gradient)(nColBins)
  exprBins <- seq(-1, 1, length=nColBins+1)
  
  nullPlot(c(0, ncol(gene.expr)), c(0, nrow(gene.expr)), xlab, ylab, cex.lab)
  
  for (jj in 1:ncol(gene.expr)) {
    for (ii in 1:nrow(gene.expr)) {
      rect(
        xleft = jj - 1,
        ybottom = nrow(gene.expr) - (ii-1),
        xright = jj,
        ytop = nrow(gene.expr) - ii,
        col=findColInGrad(gene.expr[ii, jj], exprBins, colGrad),
        border=NA
      )
    }
  }
  box()
  addTitle(main, cex.title)
}

#' Coexpression plot
#' 
#' Plot a triangle of the coexpression
#' 
#' @param coexpresion a square matrix containing the pairwise gene coexpression.
#' @param module.labels an optional named vector assigning each node to a module.
#' @param heatmap.gradient A vector of colors to use for (or interpolate over)
#'        to render the coexpression.
#' @param cex.title cex for the title text
#' @param cex.lab cex for the axis label
#' @param main title for the plot
#' @param adjacency logical; are we plotting the adjacency instead? Affects 
#'  automatic color palette and title selection
#' @param xlab label for the x axis
#'        
#' @export
plotCoexpression <- function(
  coexpression, module.labels=NULL, heatmap.gradient, adjacency=FALSE, 
  cex.title=1.4, cex.lab=1,
  main=ifelse(adjacency, "Adjacencies", "Coexpression"), xlab="genes"
) {
  coexpression <- dynamicMatLoad(coexpression)
  if (nrow(coexpression) != ncol(coexpression))
    stop("expecting a square edge.matrix!")
  
  nGenes <- nrow(coexpression)
  
  if(missing(heatmap.gradient)) {
    heatmap.gradient <- custom.palette(!adjacency)
  }
  nColBins <- 255 
  colGrad <- colorRampPalette(heatmap.gradient)(nColBins)
  
  if (adjacency) {
    edgeBins <- seq(0, 1, length=nColBins+1)
  } else {
    edgeBins <- seq(-1, 1, length=nColBins+1)
  }
  
  # create empty plot window
  nullPlot(c(0, nGenes), c(0, nGenes/2), xlab, cex.lab=cex.lab)
  
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
  addTitle(main, cex.title)
}

#-----------------------------------------------------------------
# Other useful functions not used in the current preservationPlot
#-----------------------------------------------------------------

#' Plot Module legend
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
  nullPlot(c(0, nGenes), c(0,1))
  
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
  nullPlot(c(0, nGenes), c(0, 0))
  
  axis(at=(1:nGenes)-0.5, side=1, labels=gene.labels, las=2, cex.axis=cex.axis)
  par(mar=old.mar)
}

#----------------------------
# Helper Functions
#----------------------------

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

# Because barplot has weird margins
myBarPlot <- function(
  heights, height.lim, cols=NA, bar.width=0.8, horiz=FALSE, border="black",
  cex.axis=1, cex.title, main, bar.lab, cex.lab
) {
  bars.lim <- c(0, length(heights))
  if (missing(height.lim)) {
    height.lim <- range(heights)
  }
  
  if (horiz) {
    xlim <- height.lim
    ylim <- bars.lim
  } else {
    xlim <- bars.lim
    ylim <- height.lim
  }
  
  xlim <- getLimsForManualUsr(xlim)
  ylim <- getLimsForManualUsr(ylim)
  plot(
    0, type='n', xaxt='n', yaxt='n', xlim=xlim, ylim=ylim, xlab="", ylab="", 
    frame.plot = FALSE
  )
  
  for (ii in seq_along(heights)) {
    bar.min <- max(0, height.lim[1])
    bar.offset <- (1 - bar.width)/2
    
    if (horiz) {
      xleft <- bar.min
      xright <- heights[ii]
      ybottom <- length(heights) - (ii - bar.offset)
      ytop <- length(heights) - (ii - 1 + bar.offset)
    } else {
      xleft <- ii - 1 + bar.offset
      xright <- ii - bar.offset
      ybottom <- bar.min
      ytop <- heights[ii]
    }
    if (length(cols) == 1) {
      col <- cols
    } else {
      col <- cols[ii]
    }
    rect(xleft, ybottom, xright, ytop, col=col, border=border)
  }
  
  if (horiz) {
    abline(v=0)
    # only print every third tick
    ticks <- axTicks(side=1)
    ticks <- ticks[seq(1, length(ticks), 3)]
    axis(side=1, at=ticks, labels=format(ticks, digits=2), cex.axis=cex.axis)
    mtext(bar.lab, side=2, cex=cex.lab)
  } else {
    abline(h=0)
    # only print every second tick because las=2 results in overlapping labels
    ticks <- axTicks(side=2)
    ticks = ticks[seq(1, length(ticks), 2)]
    axis(
      side=2, at=ticks, labels=format(ticks, digits=2), cex.axis=cex.axis, 
      las=2
    )
    mtext(bar.lab, side=1, cex=cex.lab)
  }
  addTitle(main, cex.title)
} 

# Create an empty plot
nullPlot <- function(xlim=c(0,0), ylim=c(0,0), xlab="", ylab="", cex.lab=1) {
  xlim <- getLimsForManualUsr(xlim)
  ylim <- getLimsForManualUsr(ylim)
  plot(
    0, type='n', xaxt='n', yaxt='n', xlim=xlim, ylim=ylim, xlab="", ylab="", 
    frame.plot = FALSE, cex.lab=cex.lab
  )
  mtext(xlab, side=1, cex=cex.lab)
  mtext(ylab, side=2, cex=cex.lab)
}

# Modified Brewer Palette
# 
# RColorBrewer palette "RdYlBu" with the middle color replaced with white. This
# gives a nicer contrast than the "RdBu" palette
# 
custom.palette <- function(diverging=TRUE) {
  cols <- rev(brewer.pal(11, "RdYlBu"))
  cols[6] <- "#FFFFFF"
  
  if (diverging) {
    cols
  } else {
    tail(cols, 6)
  }
}

# isRelative: should the expression be considered diverging around 0?
expression.palette <- function(diverging) {
  if (diverging) {
    brewer.pal(6, "PRGn")
  } else {
    rev(brewer.pal(9, "Greens"))
  }
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

addTitle <- function(title, cex) {
  mtext(title, side=3, cex=cex, line=0.4)
}
