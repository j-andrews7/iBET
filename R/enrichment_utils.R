
#' @title Dot Plot Enrichment Visualization
#' @description Creates interactive dot plot for enrichment results
#' @param enrich Enrichment result object (data.frame-like)
#' @param num.sets Number of top results to show (default: 10)
#' @param colour.by Column name for color gradient (default: "p.adjust")
#' @param title.for.plot Plot title
#' @param x.axis X-axis variable ("GeneRatio" or "Count")
#' @param text.size Y-axis text size
#' @param colour1 Low end color
#' @param colour2 High end color
#' @return plotly object
#' @importFrom rlang sym .data
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous theme_minimal element_text labs scale_color_gradient reorder
#' @importFrom plotly ggplotly
#' @author Jacob Martin
#' @export
make_dot_plot <- function(enrich, num.sets = 10, colour.by = "p.adjust", title.for.plot = "Title", 
                        x.axis = "GeneRatio", text.size = 10, colour1 = "red", colour2 = "blue", 
                        download.type = "png", font.type = "Arial"){
    enrich <- enrich[order(enrich$p.adjust), ]
    enrich <- enrich[seq(num.sets), ]
    dotPlot <- ggplot(enrich, 
        aes(x = .data[[x.axis]], y = reorder(Description, .data[[x.axis]]), color = .data[[colour.by]], size = Count, text = Description))+
        geom_point(alpha = 0.7)+
        scale_size_continuous(range = c(3, 7), name = "Count")+ 
        theme_minimal() +
        theme(
            axis.text.y = element_text(color = "black", size = text.size, angle = 0, hjust = .5, vjust = .5, face = "plain"),
            axis.text.x = element_text(color = "black", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain")
        )+
        labs(x = paste(x.axis, "(log10)"), y = "Description", title = title.for.plot)+
        scale_color_gradient(low = colour1, high = colour2, name = colour.by)
    
    plotlyOut <- ggplotly(dotPlot, tooltip = c("text", "x", "size", colour.by)) |>
        layout(title = list(text = "Click To Edit Title", font = font.type, x = 0.5, xanchor = "center", y = 0.98, yanchor = "top"),  # Editing the Plotly axis titles and download button 
               xaxis = list(title = list(text = x.axis, font = list(size = 18, familly = font.type, color = "black")))) |>
                    config(editable = TRUE, edits = list(titleText = TRUE, axisTitleText = TRUE), 
                        toImageButtonOptions = list(format = download.type, filename = "dot_plot", height = 500, width = 700, scale = 1),
                        displaylogo = FALSE) # Hiding Plotly Logo

    return(plotlyOut)

}
#' @title Bar Plot Enrichment Visualization
#' @description Creates interactive horizontal bar plot for enrichment analysis results
#' @param enrich Enrichment result object (e.g., from clusterProfiler)
#' @param num.sets Number of top results to display (default: 10)
#' @param colour.by Column name for fill color gradient (default: "p.adjust")
#' @param title.for.plot Plot title (default: "Title")
#' @param x.axis X-axis value column ("GeneRatio" or "Count", default: "GeneRatio")
#' @param text.size Axis text size (default: 10)
#' @param colour1 Low-end color for gradient (default: "red")
#' @param colour2 High-end color for gradient (default: "blue")
#' @return Interactive plotly object
#' @importFrom rlang sym .data
#' @importFrom ggplot2 ggplot aes geom_col reorder scale_fill_gradient labs theme_minimal
#' @importFrom ggplot2 element_text coord_flip
#' @importFrom plotly ggplotly
#' @author Jacob Martin
#' @export 
make_bar_plot <- function(enrich, num.sets = 10, colour.by = "p.adjust", title.for.plot = "Title", 
                        x.axis = "GeneRatio", text.size = 10, colour1 = "red", colour2 = "blue",
                        download.type = "png", font.type = "Arial"){
    
    enrich_df <- as.data.frame(enrich)[order(as.data.frame(enrich)[[colour.by]]), ][seq(num.sets), ]
    enrich_df[[x.axis]] <- as.numeric(as.character(enrich_df[[x.axis]]))
    barPlot <- ggplot(enrich_df, aes(x = reorder(Description, !!sym(x.axis)), y = !!sym(x.axis), fill = .data[[colour.by]], text = Description)) +
        geom_col(alpha = 0.8, width = 0.7) +
        scale_fill_gradient(low = colour1, high = colour2, name = colour.by) +
        labs(x = "Description", y = paste(x.axis, "(log10)"), title = title.for.plot) +
        theme_minimal() +
        theme(axis.text.x = element_text(color = "black", size = text.size, angle = 45, hjust = 1),
              axis.text.y = element_text(color = "black", size = text.size)) +
        coord_flip()  # Horizontal bars for long descriptions
    
    plotlyOut <- ggplotly(barPlot, tooltip = c("text", "y", "fill")) |> 
        layout(title = list(text = "Click To Edit Title", font = font.type, x = 0.5, xanchor = "center", y = 0.98, yanchor = "top"), # Editing the Plotly axis titles and download button 
               xaxis = list(title = list(text = x.axis, font = list(size = 18, familly = font.type, color = "black")))) |>
                    config(editable = TRUE, edits = list(titleText = TRUE, axisTitleText = TRUE),
                    toImageButtonOptions = list(format = download.type, filename = "bar_plot", height = 500, width = 700, scale = 1), 
                    displaylogo = FALSE) # Hiding Plotly Logo 
    

    return(plotlyOut)
}
#' @title Heatmap Enrichment Visualization
#' @description Creates interactive binary heatmap for enrichment analysis results showing gene presence/absence
#' @param enrich Enrichment result object (e.g., from clusterProfiler)
#' @param num.sets Number of top results to display (default: 10)
#' @param colour Tile color for gene presence (default: "black")
#' @param size.text Axis text size (default: 10)
#' @param title Title of the plot
#' @return Interactive plotly object
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual theme_minimal element_text labs
#' @importFrom plotly ggplotly
#' @author Jacob Martin
#' @export
make_binary_heatmap <- function(enrich, num.sets = 10, colour = "black", size.text = 10, title = NULL, download.type = "png", font.type = "Arial"){

    top_terms <- head(enrich, num.sets)
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
#Above shows the data.frame generated to plot the gg heat map in a binary fashion 
    plot_data <- merge(full_grid, long_df, by = c("gene", "pathway"), all.x = TRUE)
    plot_data$presence[is.na(plot_data$presence)] <- 0L
    plot_data$presence <- factor(plot_data$presence, levels = c(0, 1))

    heat <- ggplot(plot_data, aes(x = gene, y = pathway, fill = presence)) +
        geom_tile(color = "white", linewidth = 1) +
        scale_fill_manual(values = c("0" = "white", "1" = colour)) +
        theme_minimal(base_size = 11) +
        theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
        axis.text.y = element_text(size = size.text)
        ) +
        labs(x = "Genes", y = "Enriched GO Terms", fill = "Gene\nPresence", 
            title = title)
            
    heatPlotlyOut <- ggplotly(heat, tooltip = c("x", "y", "fill")) |>
        layout(title = list(text = "Click To Edit Title", font = font.type, x = 0.5, xanchor = "center", y = 0.98, yanchor = "top"),  # Editing the Plotly axis titles and download button 
            xaxis = list(title = list(text = "Genes", font = list(size = 18, familly = font.type, color = "black")))) |>
                    config(editable = TRUE, edits = list(titleText = TRUE, axisTitleText = TRUE),
                    toImageButtonOptions = list(format = download.type, filename = "heatMap_plot", height = 500, width = 700, scale = 1),
                    displaylogo = FALSE) #Hiding Plotly Logo 
    return(heatPlotlyOut)
    }
