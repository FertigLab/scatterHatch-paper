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
filePath <- "C:\\umd\\scatterHatch\\figuresPaper\\methods\\onlyColor"
plt <- ggplot2::ggplot(data=umapCoords, ggplot2::aes(x=x, y=y))
xRange <- ggplot2::ggplot_build(plt)$layout$panel_params[[1]]$x.range
yRange <- ggplot2::ggplot_build(plt)$layout$panel_params[[1]]$y.range
##dev.new(width=w, height=h, noRStudioGD = TRUE)
##lineWidth <- 0.15
##patternList <- list(list(pattern="\\", lineWidth=lineWidth*1.5, density = 0.6))
# in order of cell type 1,2,3,4
colorPalette <- c("#009E73", "#E69F00", "#F0E442", "#56B4E9")
for (i in seq(4)){
  cellTypeCoords = umapCoords[umapCoords$cluster == i, ]
  plt = ggplot(data=cellTypeCoords, aes(x=x, y=y, color=cluster)) + geom_point(size=0.6, alpha=0.5) +
    scale_color_manual(values = colorPalette[i]) + theme_void() + labs(y= "", x = "") + theme(legend.position = "none") +
    xlim(xRange) + ylim(yRange)
  dev.new(width = w, height = h, noRStudioGD = TRUE)
  plot(plt)
  ggsave(paste0(filePath, "\\cellType", as.character(i), ".png"), width = w, height = h, dpi = 1200, bg="white")
  dev.off()
}
##plt <- scatterHatch(umapCoords, "x", "y", color_by = "cluster", legendTitle = "Cluster", pointSize=1, patternList = patternList, colorPalette = colorPalette)
##plt <- plt + theme_void() + labs(y= "", x = "") + theme(legend.position = "none")