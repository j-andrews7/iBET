# library(shiny)
# library(shinydashboard)
# library(plotly)
# library(shinyBS)
# library(shinyjs)
# library(dashboardthemes)
# library(bslib)
# library(shinyjqui)
# library(clusterProfiler)
# library(ggtree)
# library(enrichplot)
# library(org.Hs.eg.db)
# library(GOSemSim)
# library(forcats)
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


#' @title Bar Plot Enrichment Visualization
#' @description Creates interactive horizontal bar plot for enrichment analysis results
#' @param enrich Enrichment result object (e.g., from clusterProfiler)
#' @param numSets Number of top results to display (default: 10)
#' @param colourBy Column name for fill color gradient (default: "p.adjust")
#' @param titleForPlot Plot title (default: "Title")
#' @param x_axis X-axis value column ("GeneRatio" or "Count", default: "GeneRatio")
#' @param textSize Axis text size (default: 10)
#' @param colour1 Low-end color for gradient (default: "red")
#' @param colour2 High-end color for gradient (default: "blue")
#' @return Interactive plotly object
#' @importFrom rlang sym .data
#' @importFrom ggplot2 ggplot aes geom_col reorder scale_fill_gradient labs theme_minimal
#' @importFrom ggplot2 element_text coord_flip
#' @importFrom plotly ggplotly
#' @author Jacob Martin
#' @export

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
#' @title Heatmap Enrichment Visualization
#' @description Creates interactive binary heatmap for enrichment analysis results showing gene presence/absence
#' @param enrich Enrichment result object (e.g., from clusterProfiler)
#' @param numSets Number of top results to display (default: 10)
#' @param colour Tile color for gene presence (default: "black")
#' @param sizeText Axis text size (default: 10)
#' @param title Title of the plot
#' @return Interactive plotly object
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual theme_minimal element_text labs
#' @importFrom plotly ggplotly
#' @author Jacob Martin
#' @export

heatPlotEnrichment <- function(enrich, numSets = 10, colour = "black", sizeText = 10, title = NULL){

  top_terms <- head(enrich, numSets)
  gene_list <- strsplit(top_terms$geneID, "/")
  names(gene_list) <- top_terms$Description
  
  long_df <- data.frame(
    pathway = rep(names(gene_list), lengths(gene_list)),
    gene    = unlist(gene_list),
    presence = 1L
  )
  
  # Rest of your code stays the same
  all_genes    <- sort(unique(long_df$gene))
  all_pathways <- sort(unique(long_df$pathway))

  full_grid <- expand.grid(
    gene    = all_genes,
    pathway = all_pathways,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  plot_data <- merge(full_grid, long_df, by = c("gene", "pathway"), all.x = TRUE)
  plot_data$presence[is.na(plot_data$presence)] <- 0L
  plot_data$presence <- factor(plot_data$presence, levels = c(0, 1))

  heat <- ggplot(plot_data, aes(x = gene, y = pathway, fill = presence)) +
    geom_tile(color = "white", linewidth = 1) +
    scale_fill_manual(values = c("0" = "white", "1" = colour)) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
      axis.text.y = element_text(size = sizeText)
    ) +
    labs(x = "Genes", y = "Enriched GO Terms", fill = "Gene\nPresence", 
         title = title)
         
  heatPlotlyOut <- ggplotly(heat, tooltip = c("x", "y", "fill"))
  return(heatPlotlyOut)
}

# test <- readRDS("/Volumes/JM/Development/iBET/Play/ExampleEnrichment.rds")
# test1 <- as.data.frame(test)

