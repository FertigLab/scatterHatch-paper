library(scatterHatch)
library(ggplot2)

# Supplementary Figure 4

pdacFrameDir = "C:\\umd\\scatterHatch\\figuresPaper\\supplementaryFigures\\suppFig3\\"

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


w = 4
h = 3
dev.new(width = w, height = h, noRStudioGD = TRUE)
patternList = list(list(pattern = "\\", density=1/4), list(pattern = "-", density = 1/4), list(pattern = "/"), list(pattern = ""))
supplementaryFigure4 = scatterHatch(data = pdacData, x = "Xt", y = "Yt", color_by = "location", legendTitle = "Tissue Type", patternList = patternList) +
    theme_void()
plot(supplementaryFigure4)
ggsave("C:\\umd\\scatterHatch\\figuresPaper\\supplementaryFigures\\Supplementary_Figure4.png", width=w, height=h, dpi = 350, bg="White")
dev.off()