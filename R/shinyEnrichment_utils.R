library(shiny)
library(shinydashboard)
library(plotly)
library(shinyBS)
library(shinyjs)
library(dashboardthemes)
library(bslib)
library(shinyjqui)
library(clusterProfiler)
library(ggtree)
library(enrichplot)
library(org.Hs.eg.db)
library(GOSemSim)
library(forcats)
#' @title Dot Plot Enrichment Visualization
#' @description Creates interactive dot plot for enrichment results
#' @param enrich Enrichment result object (data.frame-like)
#' @param numSets Number of top results to show (default: 10)
#' @param colourBy Column name for color gradient (default: "p.adjust")
#' @param titleForPlot Plot title
#' @param x_axis X-axis variable ("GeneRatio" or "Count")
#' @param textSize Y-axis text size
#' @param colour1 Low end color
#' @param colour2 High end color
#' @return plotly object
#' @importFrom rlang sym .data
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous theme_minimal element_text labs scale_color_gradient reorder
#' @importFrom plotly ggplotly
#' @aurhtor Jacob Martin
#' @export

dotPlotEnrichment <- function(enrich, numSets = 10, colourBy = "p.adjust", titleForPlot = "Title", x_axis = "GeneRatio", textSize = 10, colour1 = "red", colour2 = "blue"){
    enrich <- enrich[order(enrich$p.adjust), ]
    enrich <- enrich[1:numSets, ]
    dotPlot <- ggplot(enrich, aes(x = !!sym(x_axis), y = reorder(Description, !!sym(x_axis)), color = .data[[colourBy]], size = Count, group = Description, text = Description)) +
        geom_point(alpha = 0.7) +
        scale_size_continuous(range = c(3, 7), name = "Count") + 
        theme_minimal() +
        theme(
            axis.text.y = element_text(color = "black", size = textSize, angle = 0, hjust = .5, vjust = .5, face = "plain"),
            axis.text.x = element_text(color = "black", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain")
        )+
        labs(x = paste(x_axis, "(log10)"), y = "Description", title = titleForPlot)+
        scale_color_gradient(low = colour1, high = colour2, name = colourBy)
    plotlyOut <- ggplotly(dotPlot, tooltip = c("x", "y", "size", colourBy))
    return(plotlyOut)
}


# #'@title BarPlot Function for enrichment Results
# #' @param enrichObj enrichGO object from cluster profiler
# #' @param numSets the number of gene sets you want to be displayed in the plot
# #' @param colourBy The colour gradient is determined by: pvalue, p.adjust or qvalue 
# #' @param titleForPlot Title for the plot
# #' @param x_axis Whether the x axis represents the GeneRatio or Count 
# #' @importFrom plotly ggplotly 
# #' @importFrom clusterProfiler barplot
# #' @importFrom ggplot2 ggplot
# #' @author Jacob Martin 
# #' @export

barPlotEnrichment <- function(enrich, numSets = 10, colourBy = "p.adjust", titleForPlot = "Title", x_axis = "GeneRatio", textSize = 10, colour1 = "red", colour2 = "blue"){
    
    enrich_df <- as.data.frame(enrich)[order(as.data.frame(enrich)[[colourBy]]), ][1:numSets, ]
    enrich_df[[x_axis]] <- as.numeric(as.character(enrich_df[[x_axis]]))
    barPlot <- ggplot(enrich_df, aes(x = reorder(Description, !!sym(x_axis)), y = !!sym(x_axis), fill = .data[[colourBy]], text = Description)) +
        geom_col(alpha = 0.8, width = 0.7) +
        scale_fill_gradient(low = colour1, high = colour2, name = colourBy) +
        labs(x = "Description", y = paste(x_axis, "(log10)"), title = titleForPlot) +
        theme_minimal() +
        theme(axis.text.x = element_text(color = "black", size = textSize, angle = 45, hjust = 1),
              axis.text.y = element_text(color = "black", size = textSize)) +
        coord_flip()  # Horizontal bars for long descriptions
    
    plotlyOut <- ggplotly(barPlot, tooltip = c("text", "y", "fill"))
    return(plotlyOut)
}

# #'@title HeatMap Plot
# #' @param enrich enrichGO object from cluster profiler
# #' @param numSets the number of gene sets you want to be displayed in the plot
# #' @importFrom plotly ggplotly 
# #' @importFrom clusterProfiler heatplot
# #' @importFrom ggplot2 ggplot
# #' @author Jacob Martin 
# #' @export
# #' 
# heatPlotEnrichment <- function(enrich, numSets = 10){
#     heatPlot <- heatplot(enrich, showCategory = numSets)
#     plotlyHeatOut <- ggplotly(heatPlot, tooltip = c("text", "Count", "GeneRatio", "Description"))
#     return(plotlyHeatOut)
# }

# #'@title Tree Plot
# #' @param enrich enrichGO object from cluster profiler
# #' @param numSets the number of gene sets you want to be displayed in the plot
# #' @param colourBy The colour gradient is determined by: pvalue, p.adjust or qvalue 
# #' @importFrom plotly ggplotly 
# #' @importFrom enrichplot treeplot
# #' @importFrom ggplot2 ggplot
# #' @author Jacob Martin 
# #' @export
# treePlotEnrichment <- function(enrich, numSets = 10, colourBy = "p.adjust"){
#     treePlot <- treeplot(enrich, showCategory = numSets, color = colourBy)
#     plotlyHeatOut <- ggplotly(treePlot)
#     return(plotlyHeatOut)
# }




# dotPlotEnrichmentCustom <- ggplot(aes(x = enrich$))
#Description: 
#pvalue 
#p.adjust 
#GeneRatio 
#Count
test <- readRDS("/Volumes/JM/Development/iBET/Play/ExampleEnrichment.rds")
test1 <- as.data.frame(test)
x <- dotPlotEnrichment(enrich = test, numSets = 10, colourBy = "p.adjust", titleForPlot = "Title", x_axis = "GeneRatio", textSize = 10, colour1 = "red", colour2 = "blue")
y <- barPlotEnrichment(enrich = test, numSets = 10, colourBy = "p.adjust", titleForPlot = "Title", x_axis = "GeneRatio", textSize = 10, colour1 = "red", colour2 = "blue")

