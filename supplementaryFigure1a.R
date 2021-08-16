library(scatterHatch)
library(ggplot2)
library(viridis)
library(ggthemes)
library(scales)
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

# reordering viridis color palette so frame colors are more visually distinct
numberOfFrames = length(unique(pdacData$frame))
viridisColorPalette = viridis_pal()(numberOfFrames)
reorder = c()
for (i in seq(1, 9)){
    reorder = c(reorder, seq(i, numberOfFrames, by = 9))
}
viridisColorPalette = viridisColorPalette[reorder]

# dittoSeq color palette
dittoColors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
                "#007756", "#D5C711", "#005685", "#A04700", "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
                "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C", "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E", "#38B7FF", "#FF9B4D",
                "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45", "#AA9F0D", "#00446B", "#803800", "#8D3666", "#3D3D3D")
colorPalette = rep(dittoColors, times=ceiling(numberOfFrames/40))[1:numberOfFrames]

# Creating Figure 1A
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
ggplotColors = gg_color_hue(40)

reorder = c()
for (i in seq(1, 8)){
    reorder = c(reorder, seq(i, length(dittoColors), by = 8))
}
ggplotColors = ggplotColors[reorder]
ylabelTheme = theme(axis.title.y=element_text(size=12, angle = 90, margin=margin(0,-20,0,0)))

dev.new(width = 8, height = 8, noRStudioGD = TRUE)

# normal vision
ggplotColorsDisplay = cvdPlot(displayColors(ggplotColors, aes(x = rep(seq(1, 10), 4), y = rep(seq(1, 4), each = 10))), layout = c('origin')) + scale_y_continuous('ggplot Colors') + ylabelTheme
ggplotColorsDisplay$layers[[2]] <- NULL

dittoColorsDisplay = cvdPlot(displayColors(dittoColors, aes(x = rep(seq(1, 10), 4), y = rep(seq(1, 4), each = 10))), layout = c('origin')) + scale_y_continuous('dittoSeq Colors') + ylabelTheme
dittoColorsDisplay$layers[[2]] <- NULL
originalColors = ggarrange(ggplotColorsDisplay, dittoColorsDisplay, nrow = 2)
originalColors = annotate_figure(originalColors, top = text_grob("A", color = "black", face = "bold", size = 14, hjust = 12))

# deuteranope
ggplotColorsDisplay = cvdPlot(displayColors(ggplotColors, aes(x = rep(seq(1, 10), 4), y = rep(seq(1, 4), each = 10))), layout = c('deuteranope')) + scale_y_continuous('ggplot Colors') + ylabelTheme
ggplotColorsDisplay$layers[[2]] <- NULL

dittoColorsDisplay = cvdPlot(displayColors(dittoColors, aes(x = rep(seq(1, 10), 4), y = rep(seq(1, 4), each = 10))), layout = c('deuteranope')) + scale_y_continuous('dittoSeq Colors') + ylabelTheme
dittoColorsDisplay$layers[[2]] <- NULL
deuteranopeColors = ggarrange(ggplotColorsDisplay, dittoColorsDisplay, nrow = 2)
deuteranopeColors = annotate_figure(deuteranopeColors, top = text_grob("B", color = "black", face = "bold", size = 14, hjust = 11.25))

# protanopia
ggplotColorsDisplay = cvdPlot(displayColors(ggplotColors, aes(x = rep(seq(1, 10), 4), y = rep(seq(1, 4), each = 10))), layout = c('protanope')) + scale_y_continuous('ggplot Colors') + ylabelTheme
ggplotColorsDisplay$layers[[2]] <- NULL

dittoColorsDisplay = cvdPlot(displayColors(dittoColors, aes(x = rep(seq(1, 10), 4), y = rep(seq(1, 4), each = 10))), layout = c('protanope')) + scale_y_continuous('dittoSeq Colors') + ylabelTheme
dittoColorsDisplay$layers[[2]] <- NULL
protanopeColors = ggarrange(ggplotColorsDisplay, dittoColorsDisplay, nrow = 2)
protanopeColors = annotate_figure(protanopeColors, top = text_grob("C", color = "black", face = "bold", size = 14, hjust = 11.25))

# desaturate
ggplotColorsDisplay = cvdPlot(displayColors(ggplotColors, aes(x = rep(seq(1, 10), 4), y = rep(seq(1, 4), each = 10))), layout = c('desaturate')) + scale_y_continuous('ggplot Colors') + ylabelTheme
ggplotColorsDisplay$layers[[2]] <- NULL

dittoColorsDisplay = cvdPlot(displayColors(dittoColors, aes(x = rep(seq(1, 10), 4), y = rep(seq(1, 4), each = 10))), layout = c('desaturate')) + scale_y_continuous('dittoSeq Colors') + ylabelTheme
dittoColorsDisplay$layers[[2]] <- NULL
desaturateColors = ggarrange(ggplotColorsDisplay, dittoColorsDisplay, nrow = 2)
desaturateColors = annotate_figure(desaturateColors, top = text_grob("D", color = "black", face = "bold", size = 14, hjust = 11.25))

ggarrange(originalColors, deuteranopeColors, protanopeColors,  desaturateColors, nrow=2, ncol = 2)

ggsave("C:\\umd\\scatterHatch\\figuresPaper\\supplementaryFigures\\Supplementary_Figure1a.png", width = 8, height = 8, dpi = 96, bg = "white")
dev.off()
