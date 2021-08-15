library(umap)
library(tiff)
library(ggplot2)
library(scatterHatch)
library(colorBlindness)
library(ggpubr)

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

w <- 3
h <- 3
filePath <- "C:\\umd\\scatterHatch\\figuresPaper\\umap\\"
dev.new(width=w, height=h, noRStudioGD = TRUE)
lineWidth <- 0.15
patternList <- list(list(pattern="\\", lineWidth=lineWidth*1.5), list(pattern="|", lineWidth=lineWidth), list(pattern="-", lineWidth=lineWidth, density=1/4), list(pattern="x", lineWidth=lineWidth))
plt <- scatterHatch(umapCoords, "x", "y", factor = "cluster", legendTitle = "Cluster", pointSize=1, patternList = patternList)
plt <- plt + theme_void() + theme(legend.position = "none") + labs(y= "", x = "")
titles <- c("(A) Normal Vision", "(B) Deuteranopia", "   (C) Protanopia", "(D) Monochromacy")
perceptions <- c("origin", "deuteranope", "protanope", "desaturate")
i <- 1
for (i in seq(4)){
    colorPerception <- cvdPlot(plt, layout=perceptions[i])
    colorPerception$layers[[2]]$data$text <- titles[i]
    plot(colorPerception)
    ggsave(paste0(filePath, perceptions[i], ".tiff"), width=w, height=h, dpi=1200)
    dev.off()
    dev.new(width=w, height=h, noRStudioGD = TRUE)
}

origin <- readTIFF(paste0(filePath, "origin.tiff"))
deu <- readTIFF(paste0(filePath, "deuteranope.tiff"))
pro <- readTIFF(paste0(filePath, "protanope.tiff"))
des <- readTIFF(paste0(filePath, "desaturate.tiff"))

par(mar=rep(0,4))
outline = matrix(c(rep(c(1), times = 1), rep(2, times = 1), rep(3, times = 1), rep(4, times = 1)), 
                 ncol = 2, nrow = 2, byrow = TRUE)

layout(outline)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(origin,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(deu,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(pro,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(des,0,0,1,1)
dev.print(tiff, paste0(filePath, "umap1.tiff"), width=3, height=3, units="in", res=1200)
dev.print(tiff, paste0(filePath, "umapLosRes.tiff"), width=3, height=3, units="in", res=350)
