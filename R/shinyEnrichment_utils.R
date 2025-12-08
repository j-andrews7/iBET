
#' @title DotPlot Function For Enrichment Results 
#' @param enrich ClusterProfiler enrichment Object from enrichGO() function 
#' @param numSets the number of gene sets you want to be displayed in the plot 
#' @param colourBy Either p.adjust, qvalue or pvalue. Determines how the dots are coloured on a colour scale.
#' @param titleForPlot Title for the plot
#' @param x_axis Whether the x axis represents the GeneRatio or Count  
#' @importFrom plotly ggplotly
#' @importFrom clusterProfiler dotplot
#' @importFrom ggplot2 ggplot
#' @author Jacob Martin 
#' @export
# Dot Plot customisation: 
# Colour of the balls: p.adjust, pvalue or qvalue: 
#showCategory 
#Function for DotPlot graphing: 
# Title 
dotPlotEnrichment <- function(enrich, numSets = 10, colourBy = "p.adjust", titleForPlot = "Title", x_axis = "GeneRatio"){
    dotPlot <- dotplot(enrich, x = x_axis, showCategory = numSets, color = colourBy, title = titleForPlot)
    plotlyOut <- ggplotly(dotPlot, tooltip = c("text", "Count", "Description", colourBy))
    return(plotlyOut)
}


#'@title BarPlot Function for enrichment Results
#' @param enrichObj enrichGO object from cluster profiler
#' @param numSets the number of gene sets you want to be displayed in the plot
#' @param colourBy The colour gradient is determined by: pvalue, p.adjust or qvalue 
#' @param titleForPlot Title for the plot
#' @param x_axis Whether the x axis represents the GeneRatio or Count 
#' @importFrom plotly ggplotly 
#' @importFrom clusterProfiler barplot
#' @importFrom ggplot2 ggplot
#' @author Jacob Martin 
#' @export

barPlotEnrichment <- function(enrichObj, numSets = 10, colourBy = "p.adjust", titleForPlot = "Title", x_axis = "GeneRatio"){
    barPlot <- barplot(enrichObj, x = x_axis, showCategory = numSets, color = colourBy, title = titleForPlot, )
    plotlyBarOut <- ggplotly(barPlot, tooltip = c("text", "Count", "GeneRatio", "Description", colourBy))
    return(plotlyBarOut)
}

#'@title HeatMap Plot
#' @param enrich enrichGO object from cluster profiler
#' @param numSets the number of gene sets you want to be displayed in the plot
#' @importFrom plotly ggplotly 
#' @importFrom clusterProfiler heatplot
#' @importFrom ggplot2 ggplot
#' @author Jacob Martin 
#' @export
#' 
heatPlotEnrichment <- function(enrich, numSets = 10){
    heatPlot <- heatplot(enrich, showCategory = numSets)
    plotlyHeatOut <- ggplotly(heatPlot, tooltip = c("text", "Count", "GeneRatio", "Description"))
    return(plotlyHeatOut)
}

#'@title Tree Plot
#' @param enrich enrichGO object from cluster profiler
#' @param numSets the number of gene sets you want to be displayed in the plot
#' @param colourBy The colour gradient is determined by: pvalue, p.adjust or qvalue 
#' @importFrom plotly ggplotly 
#' @importFrom enrichplot treeplot
#' @importFrom ggplot2 ggplot
#' @author Jacob Martin 
#' @export
library(enrichplot)
treePlotEnrichment <- function(enrich, numSets = 10, colourBy = "p.adjust"){
    treePlot <- treeplot(enrich, showCategory = numSets, color = colourBy)
    plotlyHeatOut <- ggplotly(treePlot)
    return(plotlyHeatOut)
}
# test <- readRDS("/Volumes/JM/Development/iBET/Play/ExampleEnrichment.rds")
# # Before calling treePlotEnrichment()
# # test <- simplify(test, cutoff = 0.7, by = "p.adjust", select_fun = min)  # Populates @termsim
# treePlotEnrichment(test, numSets = 10)


