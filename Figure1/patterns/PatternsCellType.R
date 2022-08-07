library(umap)
library(tiff)
library(ggplot2)
library(scatterHatch)
library(colorBlindness)
library(ggpubr)

set.seed(156)

pdacData <- read.csv("C://umd//scatterHatch//scatterHatch-paper//rawdata_Figure7&8_PDAC.csv")
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

## Creating UMAPs for each cell type
w <- 3
h <- 3
filePath <- "C:\\umd\\scatterHatch\\figuresPaper\\methods\\patterns"
colorPalette <- c("#009E73", "#E69F00", "#F0E442", "#56B4E9")
lineWidth <- 0.15
patternList <- list(list(pattern="/", lineWidth=lineWidth*1.5, density = 0.6), list(pattern="x", lineWidth=lineWidth, density=0.6), list(pattern="-", lineWidth=lineWidth, density=0.5), list(pattern=""))
dev.new(width=w, height=h, noRStudioGD = TRUE)
for (i in seq(4)){
  cellTypeCoords = umapCoords[umapCoords$cluster == i, ]
  ##levels(cellTypeCoords$cluster) <- droplevels(cellTypeCoords$cluster)
  #xSpan <- range(cellTypeCoords$x); ySpan <- range(cellTypeCoords$y)
  ##gridSize = min(xSpan,ySpan)/200
  
  showOnlyCellTypeColorPalette = colorPalette
  showOnlyCellTypePatternList = patternList
  for (j in seq(4)){
    if (i != j){
      showOnlyCellTypeColorPalette[j] = "#FFFFFF"
      showOnlyCellTypePatternList[[j]]$pattern = ""
    }
  }
  
  plt = scatterHatch(data = umapCoords, x = "x", y = "y", color_by = "cluster", pointSize = 1, 
                     colorPalette = showOnlyCellTypeColorPalette, patternList = showOnlyCellTypePatternList) +  theme_void() + 
                     labs(y= "", x = "") + theme(legend.position = "none")
  
  
  plot(plt)
  ggsave(paste0(filePath, "\\cellType", as.character(i), ".png"), width = w, height = h, dpi = 1200, bg="white")
  dev.new(width=w, height=h, noRStudioGD = TRUE)
}
