library(ggplot2)
library(viridis)
library(ggthemes)
library(scales)
library(colorBlindness)
library(ggpubr)
library(scatterHatch)

pdacFrameDir = "C:\\umd\\scatterHatch\\figuresPaper\\supplementaryFigures\\suppFig2\\"

pdacData = scatterHatch::pdacData
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

# dittoSeq color palette
dittoColors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
                "#007756", "#D5C711", "#005685", "#A04700", "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
                "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C", "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E", "#38B7FF", "#FF9B4D",
                "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45", "#AA9F0D", "#00446B", "#803800", "#8D3666", "#3D3D3D")
numberOfFrames = length(unique(pdacData$frame))
colorPalette = rep(dittoColors, times=ceiling(numberOfFrames/40))[1:numberOfFrames]

# changing pattern ordering

frames = unique(pdacData$frame)
patterns = c("", "\\", "x","+", "|", "/", "x", "+", "\\", "-","|", "x", "/")
reversePatterns = rev(patterns)
reversePatterns[7] = "-"
completePatterns = c()

for (j in frames){
    row = ceiling(j/13)
    col = j - ((row-1)*13)
    
    # alternate between reverse/normal patterning every row
    if ((row %% 2) == 1){ completePatterns = c(completePatterns, patterns[col])}
    if ((row %% 2) == 0){ completePatterns = c(completePatterns, reversePatterns[col])}
}
patterns = completePatterns
patternList = vector(mode = "list", length = numberOfFrames) # initializing patternList
patternList = lapply(1:numberOfFrames, function(i){
    lineType = "solid"
    patternList[[i]] = list(pattern = patterns[i], lineWidth = 0.15, lineType = lineType)
})

# Creating PDAC by frame
w = 4
h = 3
dpi = 350

dev.new(width=w, height=h, noRStudioGD = TRUE)
plt <- scatterHatch(pdacData, "Xt", "Yt", color_by = "frame", pointSize=1, patternList = patternList, colorPalette = colorPalette)
plt <- plt + theme_void() + labs(y= "", x = "") + theme(legend.position = "none")
titles <- c("A", "B", "C", "D")
perceptions <- c("origin", "deuteranope", "protanope", "desaturate")

i <- 1
for (i in seq(4)){
    colorPerception <- cvdPlot(plt, layout=perceptions[i])
    colorPerception$layers[[2]]$data$text <- titles[i]
    plot(colorPerception)
    ggsave(paste0(pdacFrameDir, perceptions[i], ".png"), width=w, height=h, dpi=dpi, bg="white")
    dev.off()
    dev.new(width=w, height=h, noRStudioGD = TRUE, bg="white")
}

# stitching together final figure
library("png")
png("C:\\umd\\scatterHatch\\figuresPaper\\supplementaryFigures\\Supplementary_Figure2.png", width = 2400, height = 1800, bg="white")

origin <- readPNG(paste(pdacFrameDir, "origin.png", sep=""))
deu <- readPNG(paste(pdacFrameDir, "deuteranope.png", sep=""))
pro <- readPNG(paste(pdacFrameDir, "protanope.png", sep=""))
des <- readPNG(paste(pdacFrameDir, "desaturate.png", sep=""))

par(mar=rep(0,4))
outline = matrix(c(rep(1, times = 1), rep(2, times = 1), rep(3, times = 1), rep(4, times = 1)), 
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
dev.off()
#dev.print(png, "C:\\umd\\scatterHatch\\figuresPaper\\supplementaryFigures\\Supplementary_Figure2.png", width = 2400, height = 1800)