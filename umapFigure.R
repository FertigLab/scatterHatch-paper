library(umap)
library(tiff)
library(ggplot2)
library(scatterHatch)
library(colorBlindness)
library(ggpubr)
set.seed(156)

pdacData <- read.csv("C:\\umd\\summer 2020\\spatial-datasets-master\\spatial-datasets-master\\data\\2018_CyCIF_PDAC\\rawdata_Figure78_PDAC\\rawdata_Figure7&8_PDAC.csv")
dacData$cellID = paste0('cell_', 1:nrow(pdacData))
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
antiBodyData.labels <- pdacData[, "location"]
antiBodyData <- umap(antiBodyData)
umapCoords <- data.frame(data=antiBodyData$layout)
clusters <- kmeans(umapCoords, 4)
umapCoords <- data.frame(umapCoords, clusters$cluster)
colnames(umapCoords) <- c("x","y","cluster")
umapCoords$cluster <- as.factor(umapCoords$cluster)

## Creating UMAPs for each color perception
w <- 3
h <- 3
filePath <- "C:\\umd\\scatterHatch\\figuresPaper\\umap\\"
dev.new(width=w, height=h, noRStudioGD = TRUE)
lineWidth <- 0.15
patternList <- list(list(pattern="\\", lineWidth=lineWidth*1.5), list(pattern="x", lineWidth=lineWidth), list(pattern="-", lineWidth=lineWidth, density=1/4), list(pattern=""))
plt <- scatterHatch(umapCoords, "x", "y", factor = "cluster", legendTitle = "Cluster", pointSize=1, patternList = patternList)
plt <- plt + theme_void() + labs(y= "", x = "") + theme(legend.position = "none")
titles <- c("A", "B", "C", "D")
perceptions <- c("origin", "deuteranope", "protanope", "desaturate")
i <- 1
for (i in seq(4)){
    colorPerception <- cvdPlot(plt, layout=perceptions[i])
    colorPerception$layers[[2]]$data$text <- titles[i]
    plot(colorPerception)
    ggsave(paste0(filePath, "umaps\\", perceptions[i], ".tiff"), width=w, height=h, dpi=1200)
    dev.off()
    dev.new(width=w, height=h, noRStudioGD = TRUE)
}

## Creating Zoom-ins
w = 2
h = 3
bounds = umapCoords[umapCoords$x < 0.7 & umapCoords$x > -1 & umapCoords$y > -1.5 & umapCoords$y < 2,]
lineWidth = 0.3
patternList <- list(list(pattern="\\", lineWidth=lineWidth), list(pattern="x", lineWidth=lineWidth), list(pattern="-", lineWidth=lineWidth, density=1/4), list(pattern=""))
colorPalette <- colorPalette[1:4]
perceptions <- c("origin", "deuteranope", "protanope", "desaturate")

dev.new(width = w, height = h, noRStudioGD = TRUE)
plt <- scatterHatch(data = bounds, x="x", y="y", factor="cluster", pointSize=3, patternList = patternList, 
                    colorPalette = colorPalette[1:4]) + theme_void() + labs(y= "", x = "")

i <- 1
for (i in seq(4)){
    colorPerception <- cvdPlot(plt, layout=perceptions[i])
    colorPerception$layers[[2]]$data$text <- ""
    plot(colorPerception)
    ggsave(paste0(filePath, "zoomIn\\", perceptions[i], ".tiff"), width=w, height=h, dpi=1200)
    dev.off()
    dev.new(width = w, height = h, noRStudioGD = TRUE)
}

## Stitching UMAPs with zoom-in figures

origin <- readTIFF(paste0(filePath, "umaps\\origin.tiff"))
originZoom <- readTIFF(paste0(filePath, "zoomIn\\origin.tiff"))
deu <- readTIFF(paste0(filePath, "umaps\\deuteranope.tiff"))
deuZoom <- readTIFF(paste0(filePath, "zoomIn\\deuteranope.tiff"))
pro <- readTIFF(paste0(filePath, "umaps\\protanope.tiff"))
proZoom <- readTIFF(paste0(filePath, "zoomIn\\protanope.tiff"))
des <- readTIFF(paste0(filePath, "umaps\\desaturate.tiff"))
desZoom <- readTIFF(paste0(filePath, "zoomIn\\desaturate.tiff"))

par(mar=rep(0,4))
outline = matrix(c(rep(c(1), times = 3), rep(2, times = 2), rep(3, times = 3), rep(4, times = 2),
                   rep(5, times = 3), rep(6, times = 2), rep(7, times = 3), rep(8, times = 2)), 
                 ncol = 10, nrow = 2, byrow = TRUE)

layout(outline)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(origin,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(originZoom,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(deu,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(deuZoom,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(pro,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(proZoom,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(des,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(desZoom,0,0,1,1)
dev.print(tiff, paste0(filePath, "umap.tiff"), width=5, height=3, units="in", res=1200)
dev.print(tiff, paste0(filePath, "umapLowRes.tiff"), width=5, height=3, units="in", res=350)
