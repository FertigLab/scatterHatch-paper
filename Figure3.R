library(ggplot2)
library(viridis)
library(ggthemes)
library(scales)
library(colorBlindness)
library(ggpubr)
library(scatterHatch)
library(ggpubfigs)


pdacFrameDir = "C:\\umd\\scatterHatch\\figuresPaper\\supplementaryFigures\\suppFig2\\"

data(pdacData)
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
    row = ceiling(as.integer(j)/13)
    col = as.integer(j) - ((row-1)*13)
    
    # alternate between reverse/normal patterning every row
    if ((row %% 2) == 1){ completePatterns = c(completePatterns, patterns[col])}
    if ((row %% 2) == 0){ completePatterns = c(completePatterns, reversePatterns[col])}
}
patterns = completePatterns
patternList = vector(mode = "list", length = numberOfFrames) # initializing patternList
patternList = lapply(1:numberOfFrames, function(i){
    lineType = "solid"
    if (patterns[i] %in% c("+", "x")){ lineWidth = 0.08 }
    else{ lineWidth = 0.08}
    patternList[[i]] = list(pattern = patterns[i], lineWidth = lineWidth, lineType = lineType, density=1)
})

# Creating PDAC by frame

## scatterHatch
w = 4
h = 3
dpi = 750

xSpan <- abs(diff(range(pdacData$Xt))); ySpan <- abs(diff(range(pdacData$Yt)))
n <- 1.5
gridSize <- min(xSpan, ySpan)/200*n
dev.new(width=w, height=h, noRStudioGD = TRUE)

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

#colorPalette <- gg_color_hue(length(unique(pdacData$frame)))
colorPalette <- friendly_pal("muted_nine", 82, "continuous")

plt <- scatterHatch(pdacData, "Xt", "Yt", color_by = "frame", gridSize=gridSize, pointSize=0.2, patternList = patternList, colorPalette = colorPalette)
plt <- plt + ggplot2::theme_void() + ggplot2::labs(y= "", x = "") + ggplot2::theme(legend.position = "none")
titles <- c("", "", "", "")
perceptions <- c("origin", "deuteranope", "protanope", "desaturate")

i <- 1
for (i in seq(4)){
    colorPerception <- cvdPlot(plt, layout=perceptions[i])
    colorPerception$layers[[2]]$data$text <- titles[i]
    plot(colorPerception)
    ggplot2::ggsave(paste0(pdacFrameDir, perceptions[i], ".pdf"), width=w, height=h, bg="white")
    #ggplot2::ggsave(paste0(pdacFrameDir, perceptions[i], ".png"), width=w, height=h, dpi=750, bg="white")
    dev.off()
    dev.new(width=w, height=h, noRStudioGD = TRUE, bg="white")
}

## Regular Plot
w = 4
h = 3
dpi = 750
plt <- ggplot2::ggplot(data=pdacData, ggplot2::aes(x=Xt, y=Yt))
xRange <- ggplot2::ggplot_build(plt)$layout$panel_params[[1]]$x.range
yRange <- ggplot2::ggplot_build(plt)$layout$panel_params[[1]]$y.range
pdacData$frame <-  factor(pdacData$frame)
dev.new(width=w, height=h, noRStudioGD = TRUE)

plt <- ggplot2::ggplot(data=pdacData, ggplot2::aes(x=Xt, y=Yt, color=frame)) + ggplot2::geom_point(size=0.2, stroke = 0, alpha = 0.5) +
    ggplot2::theme_void() + ggplot2::labs(y= "", x = "") + ggplot2::theme(legend.position = "none") + ggplot2::xlim(xRange) + ggplot2::ylim(yRange) +
    ggplot2::scale_color_manual(values = friendly_pal("muted_nine", 82, "continuous"))

perceptions <- c("origin", "deuteranope", "protanope", "desaturate")
titles <- c("A", "B", "C", "D")

i <- 1
for (i in seq(4)){
    colorPerception <- cvdPlot(plt, layout=perceptions[i])
    colorPerception$layers[[2]]$data$text <- titles[i]
    plot(colorPerception)
    #ggplot2::ggsave(paste0(pdacFrameDir, "regular\\", perceptions[i], ".pdf"), width=w, height=h)
    ggplot2::ggsave(paste0(pdacFrameDir, "regular\\", perceptions[i], ".png"), width=w, height=h, dpi=dpi, bg="white")
    dev.off()
    dev.new(width=w, height=h, noRStudioGD = TRUE, bg="white")
}


# stitching together final figure
library("png")

origin <- readPNG(paste(pdacFrameDir, "origin.png", sep=""))
originRegular <- readPNG(paste(pdacFrameDir, "regular//origin.png", sep=""))
deu <- readPNG(paste(pdacFrameDir, "deuteranope.png", sep=""))
deuRegular <- readPNG(paste(pdacFrameDir, "regular//deuteranope.png", sep=""))
pro <- readPNG(paste(pdacFrameDir, "protanope.png", sep=""))
proRegular <- readPNG(paste(pdacFrameDir, "regular//protanope.png", sep=""))
des <- readPNG(paste(pdacFrameDir, "desaturate.png", sep=""))
desRegular <- readPNG(paste(pdacFrameDir, "regular//desaturate.png", sep=""))

par(mar=rep(0,4))
outline = matrix(c(rep(1, times = 1), rep(2, times = 1), rep(3, times = 1), rep(4, times = 1), rep(5, times = 1),
                 rep(6, times = 1), rep(7, times = 1), rep(8, times = 1)),
                 ncol = 2, nrow = 4, byrow = TRUE)

layout(outline)

plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(originRegular,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(origin,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(deuRegular,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(deu,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(proRegular,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(pro,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(desRegular,0,0,1,1)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(des,0,0,1,1)

dev.print(pdf, "C:\\umd\\scatterHatch\\figuresPaper\\supplementaryFigures\\Supplementary_Figure2.pdf", width = 6, height = 9)