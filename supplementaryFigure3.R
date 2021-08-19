library(ggplot2)
library(scatterHatch)
library(colorBlindness)
library(ggpubr)

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
numberOfFrames = length(unique(pdacData$frame))
dittoColors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
                "#007756", "#D5C711", "#005685", "#A04700", "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
                "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C", "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E", "#38B7FF", "#FF9B4D",
                "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45", "#AA9F0D", "#00446B", "#803800", "#8D3666", "#3D3D3D")
colorPalette = rep(dittoColors, times=ceiling(numberOfFrames/40))[1:numberOfFrames]

dev.new(width = 8, height = 6, noRStudioGD = TRUE)
plt = ggplot(data = pdacData, aes(x = Xt, y = Yt, color = as.factor(frame))) + geom_point(size = 1, alpha = 0.5) + scale_color_manual(values = colorPalette) + 
    theme_void() + labs(y= "", x = "") + theme(legend.position = "none")

plt = cvdPlot(plt)
plt$layers[[2]]$data$text = "A"
plt$layers[[4]]$data$text = "B"
plt$layers[[6]]$data$text = "C"
plt$layers[[8]]$data$text = "D"
plot(plt)
ggsave("C:\\umd\\scatterHatch\\figuresPaper\\supplementaryFigures\\Supplementary_Figure3.png", width = 8, height = 6, dpi = 96, bg="white")
dev.off()
