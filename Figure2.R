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

## Creating UMAPs for each color perception

### scatterHatch
w <- 3
h <- 3
n <- 1

filePath <- "C:\\umd\\scatterHatch\\figuresPaper\\umap\\"
dev.new(width=w, height=h, noRStudioGD = TRUE)
lineWidth <- 0.15
patternList <- list(list(pattern="/", lineWidth=lineWidth*1.5, density = 0.6), list(pattern="x", lineWidth=lineWidth, density=0.6), list(pattern="-", lineWidth=lineWidth, density=0.5), list(pattern=""))
colorPalette <- c("#009E73", "#E69F00", "#F0E442", "#56B4E9")

xSpan <- abs(diff(range(umapCoords$x))); ySpan <- abs(diff(range(umapCoords$y)))
gridSize <- min(xSpan,ySpan)/200

plt <- scatterHatch(umapCoords, "x", "y", color_by = "cluster", legendTitle = "Cluster", pointSize=1, patternList = patternList, colorPalette = colorPalette)
plt <- plt + theme_void() + labs(y= "", x = "") + theme(legend.position = "none")
titles <- c("", "", "", "")
perceptions <- c("origin", "deuteranope", "protanope", "desaturate")
i <- 1
for (i in seq(4)){
    colorPerception <- cvdPlot(plt, layout=perceptions[i])
    colorPerception$layers[[2]]$data$text <- titles[i]
    plot(colorPerception)
    ggsave(paste0(filePath, "umaps\\", perceptions[i], ".tiff"), width=w, height=h, dpi=1200, bg = "white")
    dev.off()
    dev.new(width=w, height=h, noRStudioGD = TRUE)
}

### regular scatterplot
w <- 3
h <- 3
i <- 1
plt <- ggplot2::ggplot(data=umapCoords, ggplot2::aes(x=x, y=y))
xRange <- ggplot2::ggplot_build(plt)$layout$panel_params[[1]]$x.range
yRange <- ggplot2::ggplot_build(plt)$layout$panel_params[[1]]$y.range
colorPalette <- c("#009E73", "#E69F00", "#F0E442", "#56B4E9")

dev.new(width=w, height=h, noRStudioGD = TRUE)
plt <- ggplot(data=umapCoords, aes(x=x, y=y, color=cluster)) + geom_point(size=1, stroke=0, alpha=0.5)
plt <- plt + theme_void() + scale_color_manual(values = colorPalette) + 
   labs(y= "", x = "") + theme(legend.position = "none") + xlim(xRange) + ylim(yRange)
titles <- c("A", "B", "C", "D")
perceptions <- c("origin", "deuteranope", "protanope", "desaturate")
for (i in seq(4)){
    colorPerception <- cvdPlot(plt, layout=perceptions[i])
    colorPerception$layers[[2]]$data$text <- titles[i]
    plot(colorPerception)
    ggsave(paste0(filePath, "umaps\\regular\\", perceptions[i], ".tiff"), width=w, height=h, dpi=1200, bg="white")
    dev.off()
    dev.new(width=w, height=h, noRStudioGD = TRUE)
}

## Creating Zoom-ins

### scatterHatch
w = 2
h = 3
bounds = umapCoords[umapCoords$x < 0.7 & umapCoords$x > -1 & umapCoords$y > -1.5 & umapCoords$y < 2,]
lineWidth = 0.3
patternList <- list(list(pattern="\\", lineWidth=lineWidth, density = 0.2), list(pattern="x", lineWidth=lineWidth, density=0.6), list(pattern="-", lineWidth=lineWidth, density=0.5), list(pattern=""))
perceptions <- c("origin", "deuteranope", "protanope", "desaturate")

dev.new(width = w, height = h, noRStudioGD = TRUE)
plt <- scatterHatch(data = bounds, x="x", y="y", color_by="cluster", pointSize=3, patternList = patternList, 
                    colorPalette = colorPalette) + theme_void() + labs(y= "", x = "")

i <- 1
for (i in seq(4)){
    colorPerception <- cvdPlot(plt, layout=perceptions[i])
    colorPerception$layers[[2]]$data$text <- ""
    plot(colorPerception)
    ggsave(paste0(filePath, "zoomIn\\", perceptions[i], ".tiff"), width=w, height=h, dpi=1200)
    dev.off()
    dev.new(width = w, height = h, noRStudioGD = TRUE)
}


### regular plot
w = 2
h = 3
bounds = umapCoords[umapCoords$x < 0.7 & umapCoords$x > -1 & umapCoords$y > -1.5 & umapCoords$y < 2,]
plt <- ggplot2::ggplot(data=bounds, ggplot2::aes(x=x, y=y))
xRange <- ggplot2::ggplot_build(plt)$layout$panel_params[[1]]$x.range
yRange <- ggplot2::ggplot_build(plt)$layout$panel_params[[1]]$y.range

colorPalette <- c("#009E73", "#E69F00", "#F0E442", "#56B4E9")

dev.new(width = w, height = h, noRStudioGD = TRUE)
plt <- ggplot(data=bounds, aes(x=x, y=y, color=cluster)) + geom_point(size=3, stroke=0, alpha=0.5) +
    scale_color_manual(values = colorPalette) + labs(y= "", x = "") + xlim(xRange) + ylim(yRange) + theme_void() +
    theme(legend.title=element_blank()) + guides(color = guide_legend(override.aes = list(size=6, alpha=1)))

i <- 1
for (i in seq(4)){
    colorPerception <- cvdPlot(plt, layout=perceptions[i])
    colorPerception$layers[[2]]$data$text <- ""
    plot(colorPerception)
    ggsave(paste0(filePath, "zoomIn\\regular\\", perceptions[i], ".tiff"), width=w, height=h, dpi=1200)
    dev.off()
    dev.new(width = w, height = h, noRStudioGD = TRUE)
}


## Stitching UMAPs with zoom-in figures

origin <- readTIFF(paste0(filePath, "umaps\\origin.tiff"))
originRegular <- readTIFF(paste0(filePath, "umaps\\regular\\origin.tiff"))
originLeg <- readTIFF(paste0(filePath, "zoomIn\\origin.tiff"))
originRegularLeg <- readTIFF(paste0(filePath, "zoomIn\\regular\\origin.tiff"))

deu <- readTIFF(paste0(filePath, "umaps\\deuteranope.tiff"))
deuRegular <- readTIFF(paste0(filePath, "umaps\\regular\\deuteranope.tiff"))
deuLeg <- readTIFF(paste0(filePath, "zoomIn\\deuteranope.tiff"))
deuRegularLeg <- readTIFF(paste0(filePath, "zoomIn\\regular\\deuteranope.tiff"))

pro <- readTIFF(paste0(filePath, "umaps\\protanope.tiff"))
proRegular <- readTIFF(paste0(filePath, "umaps\\regular\\protanope.tiff"))
proLeg <- readTIFF(paste0(filePath, "zoomIn\\protanope.tiff"))
proRegularLeg <- readTIFF(paste0(filePath, "zoomIn\\regular\\protanope.tiff"))

des <- readTIFF(paste0(filePath, "umaps\\desaturate.tiff"))
desRegular <- readTIFF(paste0(filePath, "umaps\\regular\\desaturate.tiff"))
desLeg<- readTIFF(paste0(filePath, "zoomIn\\desaturate.tiff"))
desRegularLeg <- readTIFF(paste0(filePath, "zoomIn\\regular\\desaturate.tiff"))

par(mar=rep(0,4))
outline = matrix(c(rep(1, times = 3), rep(2, times = 2), rep(3, times = 1), rep(4, times = 3), rep(5, times = 2),
                   rep(6, times = 1), rep(7, times = 3), rep(8, times = 2), rep(9, times = 1), rep(10, times = 3), rep(11, times = 2), rep(12, times = 1),
                   rep(13, times = 3), rep(14, times = 2), rep(15, times = 1), rep(16, times = 3), rep(17, times = 2),
                   rep(18, times = 1), rep(19, times = 3), rep(20, times = 2), rep(21, times = 1), rep(22, times = 3), rep(23, times = 2), rep(24, times = 1)), 
                 ncol = 12, nrow = 4, byrow = TRUE)

layout(outline)
zoomRows = 1:3600
zoomCols = 1:1850
zoomL = 0
zoomB = 0
zoomR = 1
zoomT = 1

legRows = 200:3400
legCols = 1750:2400
legL = 0.15
legB = 0
legR = 0.85
legT = 1

plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(originRegular,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(originRegularLeg[zoomRows,zoomCols,],zoomL,zoomB,zoomR,zoomT)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(originRegularLeg[legRows,legCols,],legL,legB,legR,legT)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(origin,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(originLeg[zoomRows,zoomCols,],zoomL,zoomB,zoomR,zoomT)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(originLeg[legRows,legCols,],legL,legB,legR,legT)


plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(deuRegular,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(deuRegularLeg[zoomRows,zoomCols,],zoomL,zoomB,zoomR,zoomT)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(deuRegularLeg[legRows,legCols,],legL,legB,legR,legT)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(deu,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(deuLeg[zoomRows,zoomCols,],zoomL,zoomB,zoomR,zoomT)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(deuLeg[legRows,legCols,],legL,legB,legR,legT)

plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(proRegular,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(proRegularLeg[zoomRows,zoomCols,],zoomL,zoomB,zoomR,zoomT)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(proRegularLeg[legRows,legCols,],legL,legB,legR,legT)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(pro,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(proLeg[zoomRows,zoomCols,],zoomL,zoomB,zoomR,zoomT)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(proLeg[legRows,legCols,],legL,legB,legR,legT)

plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(desRegular,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(desRegularLeg[zoomRows,zoomCols,],zoomL,zoomB,zoomR,zoomT)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(desRegularLeg[legRows,legCols,],legL,legB,legR,legT)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(des,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(desLeg[zoomRows,zoomCols,],zoomL,zoomB,zoomR,zoomT)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(desLeg[legRows,legCols,],legL,legB,legR,legT)


#dev.print(tiff, paste0(filePath, "umap.tiff"), width=5, height=6, units="in", res=1200)
dev.print(pdf, paste0(filePath, "umap.pdf"), width=5, height=6)
#dev.print(tiff, paste0(filePath, "umapLowRes.tiff"), width=5, height=6, units="in", res=350)
