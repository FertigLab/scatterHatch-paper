library(umap)
library(tiff)
library(ggplot2)
library(scatterHatch)
library(colorBlindness)
library(ggpubr)

set.seed(156)

pdacData <- read.csv("C:\\umd\\summer 2020\\spatial-datasets-master\\spatial-datasets-master\\data\\2018_CyCIF_PDAC\\rawdata_Figure78_PDAC\\rawdata_Figure7&8_PDAC.csv")
pdacData$cellID = paste0('cell_', 1:nrow(pdacData))
pdacData$Yt <- -pdacData$Yt
pancreas_frames = c(1:6, 27:31, 15:19, 40:44)
PDAC_frames = c(23:26, 35:37, 51:52, 64:65, 77)
small_intestines_frames = c(49:50, 63, 75:76, 88:89, 100:103, 112:116, 125:129, 137:140)
annotateLocation <- function(frame){
  if (frame %in% pancreas_frames){return("Pancreas")}
  if (frame %in% PDAC_frames){return("PDAC")}
  if (frame %in% small_intestines_frames){return("Small Intestine")}
  return("Other")
}
pdacData$location = sapply(pdacData$frame, annotateLocation)
random <- pdacData[sample(nrow(pdacData), 10000),]
antiBodyData <- random[,1:31]
antiBodyData.labels <- random[, "location"]
antiBodyData <- umap(antiBodyData)
umapCoords <- data.frame(data=antiBodyData$layout)
clusters <- kmeans(umapCoords, 4)
umapCoords <- data.frame(umapCoords, clusters$cluster)
colnames(umapCoords) <- c("x","y","cluster")
umapCoords$cluster <- as.factor(umapCoords$cluster)
umapCoords$location <-  antiBodyData.labels

## Necessary funcs
rotateCoords <- function(x, y, angle){
  radians <- (angle/180) * pi
  xRotate <- (x * cos(radians)) - (y * sin(radians))
  yRotate <- (x * sin(radians)) + (y * cos(radians))
  return(data.frame(x = xRotate, y = yRotate))
}

countGridPoints <- function(x, y, gridSize){
  xgridStep <- gridSize
  xRange = c(min(x)-0.5*xgridStep,max(x)+xgridStep)
  yRange = c(min(y)-0.5*gridSize,max(y)+gridSize)
  xIntervals <- as.numeric(cut(x, breaks=seq(xRange[1], xRange[2], 
                                             by=xgridStep)))
  yIntervals <- as.numeric(cut(y, breaks=seq(yRange[1], yRange[2], 
                                             by=1*gridSize)))
  
  ## ensures as row index increases, y decreases in matrix
  #yIntervals <- ceiling(diff(yRange)/gridSize)-yIntervals
  intervals <- as.data.frame(cbind(xIntervals, yIntervals))
  pointsToGrid <- as.data.frame(cbind(xIntervals, yIntervals, x, y))
  freqs <- plyr::count(intervals, vars=c("xIntervals", "yIntervals"))
  freqMat <- matrix(0, nrow=max(yIntervals), ncol=max(xIntervals))
  for (ind in seq(nrow(freqs))){
    freqMat[freqs$yIntervals[ind], freqs$xIntervals[ind]] <- freqs$freq[ind]
  }
  
  xLevels <- seq(xRange[1], xRange[2], by=xgridStep)
  yLevels <- seq(yRange[1], yRange[2], by=1*gridSize)
  return(list(xLevels, yLevels, freqMat, pointsToGrid))
}

sparsityAnnotate <- function(pointsToGrid, pointSize, xRange, yRange, whichAxis) {
  
  if (whichAxis == "x") {
    pointRadius <- abs(convertSizeToCartesian(max(pointSize,1), xRange, "x"))
  }
  if (whichAxis == "y") {
    pointRadius <- abs(convertSizeToCartesian(max(pointSize,1), yRange, "y"))
  }
  
  sparsityDistance <- pointRadius * 2  # a point away
  
  ## distance from the second closest point
  pointsToGrid$closest2Points <- spatstat.geom::nndist(pointsToGrid$x, 
                                                       pointsToGrid$y, k = 2)
  
  ## distance from the fifth closest point
  pointsToGrid$closest5Points <- spatstat.geom::nndist(pointsToGrid$x, 
                                                       pointsToGrid$y, k = 5)
  
  ## distance from the twentieth closest point
  pointsToGrid$closest20Points <- spatstat.geom::nndist(pointsToGrid$x, 
                                                        pointsToGrid$y, k = 20)
  
  ## outlying points within a group
  pointsToGrid$sparsePoints <- pointsToGrid$closest2Points > sparsityDistance
  
  ## smaller clusters within a group
  pointsToGrid$smallClusters <- (pointsToGrid$closest5Points>pointRadius*4) & 
    !pointsToGrid$sparsePoints
  return(pointsToGrid)
}

getIrregularPoints <- function(pointsToGrid, freqMat, sparsePoints, 
                               rotatedxRange, rotatedyRange, pointSize){
  
  sparsityAnnotateOutput <- sparsityAnnotate(pointsToGrid, pointSize, 
                                             rotatedxRange, rotatedyRange, 'x')
  
  ## if sparse points are not given
  if (is.null(sparsePoints)){ 
    sparsePoints <- sparsityAnnotateOutput$sparsePoints
  }
  
  sparsePointsToGrid <- pointsToGrid[sparsePoints, ]
  smallClusterToGrid <- pointsToGrid[sparsityAnnotateOutput$smallClusters, ]
  
  ## removes sparse and small cluster points from regular pattern drawing
  pointsToGrid <- pointsToGrid[!sparsePoints & 
                                 !sparsityAnnotateOutput$smallClusters, ]
  
  ## removes sparse points from 2D frequency matrix
  allIrregularPoints <- rbind(sparsePointsToGrid, smallClusterToGrid)
  for (i in seq(nrow(allIrregularPoints))){ 
    freqMat[allIrregularPoints$yIntervals[i], allIrregularPoints$xIntervals[i]] = freqMat[allIrregularPoints$yIntervals[i], allIrregularPoints$xIntervals[i]] - 1
  }
  
  return(list(sparsePointsToGrid, smallClusterToGrid, pointsToGrid, freqMat))
}

convertSizeToCartesian <- function(size, dimRange, whichAxis){
  #fontSize = size*ggplot2::.pt + ggplot2::.stroke*0.5/2
  fontSize <- size*ggplot2::.pt
  aspectRatio <- dev.size()[1]/dev.size()[2]
  if (aspectRatio >= 1){
    if (whichAxis == 'x'){
      cartesianConvert <- grid::convertWidth(
        grid::unit(fontSize, "points"), unitTo="npc", valueOnly=TRUE) * 
        diff(dimRange)/aspectRatio
    }
    if (whichAxis == 'y'){
      cartesianConvert <- grid::convertHeight(
        grid::unit(fontSize, "points"), unitTo="npc", valueOnly=TRUE) * 
        diff(dimRange)/aspectRatio
    }
  }
  
  if (aspectRatio < 1){
    if (whichAxis == 'x'){
      cartesianConvert <- grid::convertWidth(
        grid::unit(fontSize, "points"), unitTo="npc", valueOnly=TRUE) * 
        diff(dimRange)*aspectRatio
    }
    if (whichAxis == 'y'){
      cartesianConvert <- grid::convertHeight(
        grid::unit(fontSize, "points"), unitTo="npc", valueOnly=TRUE) * 
        diff(dimRange)*aspectRatio
    }
  }
  
  return(cartesianConvert/2)
}

regularPatternDraw <- function(freqMat, pointsToGrid, yBins){
  xStart <- yStart <- xEnd <- yEnd <- c()
  for (rowNum in seq(nrow(freqMat))){ # iterates by every rowNum
    rowPoints <- pointsToGrid[pointsToGrid$yIntervals == rowNum, ]
    
    yLevels <- yBins[rowNum] + diff(yBins)[1]/2 # atul's version
    if (rowNum == nrow(freqMat)){ # for bottom rowNum exception
      yLevels <- yBins[rowNum] + diff(yBins)[1]/2
    }
    
    prevCol <- 0
    lineDraw <- FALSE # whether line being drawn or not
    
    for (colNum in seq(ncol(freqMat))){
      ## starting a line segment
      if (prevCol == 0 & freqMat[rowNum, colNum]!=0){ 
        gridPoints <- rowPoints[rowPoints$xIntervals == colNum, ]
        xStart <- c(xStart, min(gridPoints$x))
        yStart <- c(yStart, yLevels)
        lineDraw <- TRUE
      }
      
      ## ending line segment (added logic for handling end of freqMat before end of points)
      if (lineDraw & (freqMat[rowNum, colNum]==0||colNum==ncol(freqMat))){ 
        if (freqMat[rowNum, colNum]==0)
          gridPoints <- rowPoints[rowPoints$xIntervals == colNum-1, ]
        else
          gridPoints <- rowPoints[rowPoints$xIntervals == colNum, ]
        xEnd <- c(xEnd, max(gridPoints$x))
        yEnd <- c(yEnd, yLevels)
        lineDraw <- FALSE
      }
      
      prevCol <- freqMat[rowNum, colNum]
    }
  }
  
  return(data.frame(xStart=xStart, xEnd=xEnd, yStart=yStart, yEnd=yEnd))
}

## Creating UMAPs for each cell type
w <- 3
h <- 3
x <- umapCoords$x
y <- umapCoords$y
filePath <- "C:\\umd\\scatterHatch\\figuresPaper\\methods\\steps\\"
colorPalette <- c("#009E73", "#E69F00", "#F0E442", "#56B4E9")
lineWidth <- 0.15
patternList <- list(list(pattern="/", lineWidth=lineWidth*1.5, density = 0.6), list(pattern="x", lineWidth=lineWidth, density=0.6), list(pattern="-", lineWidth=lineWidth, density=0.5), list(pattern=""))
dev.new(width=3, height=3, noRStudioGD = TRUE)


xSpan <- abs(diff(range(x))); ySpan <- abs(diff(range(y)))
gridSize <-  min(xSpan,ySpan)/200

plt <- ggplot2::ggplot(data=umapCoords, ggplot2::aes(x=x, y=y))
xRange <- ggplot2::ggplot_build(plt)$layout$panel_params[[1]]$x.range
yRange <- ggplot2::ggplot_build(plt)$layout$panel_params[[1]]$y.range
cellType = 3
a = 0
xGroup = umapCoords[umapCoords$cluster == cellType, ]$x
yGroup = umapCoords[umapCoords$cluster == cellType, ]$y

rotatedCoords <- rotateCoords(xGroup, yGroup, angle=a) 
rotatedCoordsRange <- rotateCoords(c(xRange[1], xRange[1], xRange[2], 
                                     xRange[2]), c(yRange[1], yRange[2],
                                                   yRange[1], yRange[2]), a)


groupGridSize <- gridSize/0.6
rotatedgridOutput <- countGridPoints(rotatedCoords$x, rotatedCoords$y, gridSize=groupGridSize)
xBins <- rotatedgridOutput[[1]]; yBins <- rotatedgridOutput[[2]]
freqMat <- rotatedgridOutput[[3]]; pointsToGrid <- rotatedgridOutput[[4]]


pointClassification <- getIrregularPoints(pointsToGrid, freqMat, sparsePoints=NULL, 
                                          xRange, yRange, pointSize=1)
sparsePointsToGrid <- pointClassification[[1]] 
smallClusterToGrid <- pointClassification[[2]]
pointsToGrid <- pointClassification[[3]]

umapCoords$type = "NotCurrentCellType"
count = 1
for (i in seq(nrow(umapCoords))){
  if(umapCoords$cluster[i] == cellType){
    umapCoords$type[i] = "Dense"
    if (count %in% as.numeric(rownames(sparsePointsToGrid))){
      umapCoords$type[i] = "Sparse"
    }
    
    if (count %in% as.numeric(rownames(smallClusterToGrid))){
      umapCoords$type[i] = "SmallCluster"
    }
    
    count = count + 1
  }
}



# plotting dense/sparse only by color

# highlight sparse/small cluster
ggplot(data=umapCoords, aes(x=x, y=y, color=type)) + geom_point(size=1, alpha=0.5, stroke=0) + 
  theme_void() + labs(y= "", x = "") + theme(legend.position = "none") +
  scale_color_manual(values = c("#666666",
                                "#FFFFFF",
                                "#F0E442",
                                "#F0E442"))
ggsave(paste0(filePath, "cellType", as.character(cellType), "ColorOnlySparseHighlight.png"), width = w, height = h, dpi = 1200, bg="white")

# highlight dense
ggplot(data=umapCoords, aes(x=x, y=y, color=type)) + geom_point(size=1, alpha=0.5, stroke=0) + 
  theme_void() + labs(y= "", x = "") + theme(legend.position = "none") +
  scale_color_manual(values = c("#F0E442",
                                "#FFFFFF",
                                "#666666",
                                "#666666"))
ggsave(paste0(filePath, "cellType", as.character(cellType), "ColorOnlyDenseHighlight.png"), width = w, height = h, dpi = 1200, bg="white")


# calculating line segments for sparse/dense

## large clusters lines
denseLines <- regularPatternDraw(freqMat, pointsToGrid, yBins)

## sparse points lines
sparseLines <- sparsePointsToGrid[,c('x','x','y','y')]
names(sparseLines) <- c("xStart","xEnd","yStart","yEnd")

## converting back to unrotated axes
## adjusting lines based on sizes of a point
adjX <- convertSizeToCartesian(1, xRange, 'x')
adjY <- convertSizeToCartesian(1, yRange, 'y')

## rotation matrix
R <- matrix(c(cos(a/180 * pi),sin(a/180 * pi),
              -sin(a/180 * pi),cos(a/180 * pi)), 2, 2) 
## rotating adjustment based on angle
rotatedAdjX <- sqrt(diag(
  R %*% diag(c(adjX, adjY), 2, 2)^2 %*% t(R)))[1]
if (a == 0){ rotatedAdjX <- adjX}
## converting back to regular coordinates
rotatedDenseStartPoints <- rotateCoords(denseLines$xStart - rotatedAdjX, 
                                   denseLines$yStart, -a)
rotatedDenseEndPoints <- rotateCoords(denseLines$xEnd + rotatedAdjX, 
                                 denseLines$yEnd, -a)
## adding line segments to plot
denseLines$xStart <- rotatedDenseStartPoints$x
denseLines$yStart <- rotatedDenseStartPoints$y
denseLines$xEnd <- rotatedDenseEndPoints$x
denseLines$yEnd <- rotatedDenseEndPoints$y




plt <- ggplot(data=umapCoords, aes(x=x, y=y, color=type)) + geom_point(size=1, alpha=0.5, stroke=0) + 
  theme_void() + labs(y= "", x = "") + theme(legend.position = "none") +
  scale_color_manual(values = c("#F0E442",
                                "#FFFFFF",
                                "#666666",
                                "#666666")) +  ggplot2::geom_segment(data=denseLines,
                                   ggplot2::aes(x=xStart, y=yStart,
                                                xend=xEnd, yend=yEnd),
                                   alpha=1,
                                   size=0.15,
                                   color="black")
plot(plt)
ggsave(paste0(filePath, "cellType", as.character(cellType), "PatternDenseHighlight.png"), width = w, height = h, dpi = 1200, bg="white")

# highlighting sparse/smallCluster points with color/pattern

## converting back to regular coordinates
rotatedSparseStartPoints <- rotateCoords(sparseLines$xStart - rotatedAdjX, 
                                        sparseLines$yStart, -a)
rotatedSparseEndPoints <- rotateCoords(sparseLines$xEnd + rotatedAdjX, 
                                      sparseLines$yEnd, -a)
## adding line segments to plot
sparseLines$xStart <- rotatedSparseStartPoints$x
sparseLines$yStart <- rotatedSparseStartPoints$y
sparseLines$xEnd <- rotatedSparseEndPoints$x
sparseLines$yEnd <- rotatedSparseEndPoints$y

## creating a unique identifier for each grid
smallClusterToGrid$gridNum <- (smallClusterToGrid$yIntervals - 1)/max(smallClusterToGrid$yIntervals) + smallClusterToGrid$xIntervals

smallClusterLines = data.frame(matrix(ncol = 4, nrow = 0))
names(smallClusterLines) <- c("xStart","xEnd","yStart","yEnd")
for (gridNum in unique(smallClusterToGrid$gridNum)){
  xRange <- smallClusterToGrid$x[smallClusterToGrid$gridNum == gridNum]
  yRange <- smallClusterToGrid$y[smallClusterToGrid$gridNum == gridNum]
  clusterEnds <- data.frame(xStart=min(xRange),xEnd=max(xRange),yStart=min(yRange),yEnd=max(yRange))
  smallClusterLines <- rbind(smallClusterLines,clusterEnds)
}

## converting back to regular coordinates
rotatedSmallClusterStartPoints <- rotateCoords(smallClusterLines$xStart - rotatedAdjX, 
                                         smallClusterLines$yStart, -a)
rotatedSmallClusterEndPoints <- rotateCoords(smallClusterLines$xEnd + rotatedAdjX, 
                                       smallClusterLines$yEnd, -a)
## adding line segments to plot
smallClusterLines$xStart <- rotatedSmallClusterStartPoints$x
smallClusterLines$yStart <- rotatedSmallClusterStartPoints$y
smallClusterLines$xEnd <- rotatedSmallClusterEndPoints$x
smallClusterLines$yEnd <- rotatedSmallClusterEndPoints$y

# adding small cluster lines to sparse lines
sparseLines$xStart <- rbind(sparseLines$xStart, smallClusterLines$xStart)
sparseLines$xEnd <- rbind(sparseLines$xEnd, smallClusterLines$xEnd)
sparseLines$yStart <- rbind(sparseLines$yStart, smallClusterLines$yStart)
sparseLines$yEnd <- rbind(sparseLines$yEnd, smallClusterLines$yEnd)

plt <- ggplot(data=umapCoords, aes(x=x, y=y, color=type)) + geom_point(size=1, alpha=0.5, stroke=0) + 
  theme_void() + labs(y= "", x = "") + theme(legend.position = "none") +
  scale_color_manual(values = c("#666666",
                                "#FFFFFF",
                                "#666666",
                                "#F0E442")) +  ggplot2::geom_segment(data=sparseLines,
                                                                     ggplot2::aes(x=xStart, y=yStart,
                                                                                  xend=xEnd, yend=yEnd),
                                                                     alpha=1,
                                                                     size=0.15,
                                                                     color="black")
plot(plt)
ggsave(paste0(filePath, "cellType", as.character(cellType), "PatternSparseHighlight.png"), width = w, height = h, dpi = 1200, bg="white")
