#' Create fgsea enrichment plot for single pathway
#'
#' @param pathway Vector of gene identifiers in pathway.
#' @param stats Named vector of gene rank statistics.
#'   Each element should be named with a gene identifier.
#' @param pathway.name Optional string to use as plot title and to access fgsea stats.
#' @param fgsea.res Optional data.frame containing fgsea results as returned by `fgsea`.
#'   If provided with \code{pathway.name}, the stats for the pathway can be included on the plot.
#' @param plot.stats Boolean indicating whether to plot the fgsea stats on the plot if
#'   \code{fgsea.res} is provided and \code{pathway.name} are provided.
#' @param dge.res Optional data.frame containing differential expression results.
#'   If provided, the differentially expressed genes will be highlighted in the rugplot.
#' @param lfc.term Column name in \code{dge.res} containing log fold change values.
#'   "auto" will attempt to automatically determine the column name.
#' @param sig.term Column name in \code{dge.res} containing significance values.
#'   "auto" will attempt to automatically determine the column name.
#' @param exp.term Column name in \code{dge.res} containing expression values.
#'   "auto" will attempt to automatically determine the column name.
#' @param id.term Column name in \code{dge.res} containing gene identifiers.
#'   "rownames" will use the rownames of \code{dge.res}.
#' @param lfc.thresh Numeric value for log fold change threshold to consider
#'   a gene differentially expressed.
#' @param sig.thresh Numeric value for significance threshold to consider
#'   a gene differentially expressed.
#' @param exp.thresh Numeric value for expression threshold to consider
#'   a gene differentially expressed.
#' @param dge.up.color Color to use for ticks in rugplot for upregulated genes.
#' @param dge.down.color Color to use for ticks in rugplot for downregulated genes.
#' @param tick.color Color to use for ticks in rugplot.
#' @param gseaParam Numeric value for GSEA parameter as used in `fgsea`.
#' @return A plotly plot.
#'
#' @importFrom fgsea plotEnrichmentData
#' @importFrom plotly ggplotly config layout add_annotations %>%
#' @importFrom ggplot2 geom_line geom_segment aes theme element_blank element_line geom_hline
#'   labs scale_color_identity geom_ribbon
#'
#' @author Jared Andrews
#' @export
plot_enrichment <- function(pathway.genes,
                            stats,
                            pathway.name = NULL,
                            fgsea.res = NULL,
                            plot.stats = TRUE,
                            dge.res = NULL,
                            lfc.term = "auto",
                            sig.term = "auto",
                            exp.term = "auto",
                            id.term = "rownames",
                            lfc.thresh = 0,
                            sig.thresh = 0.05,
                            exp.thresh = 0,
                            dge.up.color = "red",
                            dge.down.color = "blue",
                            tick.color = "black",
                            gseaParam = 1) {
    # Parameter validation
    # TODO: move this to a separate function
    if (!is.null(dge.res)) {
        dge.cols <- colnames(dge.res)

        if (lfc.term == "auto") {
            if (!any(dge.cols %in% c("log2FoldChange", "logFC", "LFC"))) {
                stop("Cannot determine significance term, please provide the column name to lfc.term")
            } else {
                lfc.term <- dge.cols[dge.cols %in% c("log2FoldChange", "logFC", "LFC")]
                # If multiple matches, just use first
                if (length(lfc.term) > 1) {
                    lfc.term <- lfc.term[1]
                }
            }
        }

        if (sig.term == "auto") {
            if (!any(dge.cols %in% c("padj", "FDR", "svalue", "adj.P.Val"))) {
                stop("Cannot determine significance term, please provide the column name to sig.term")
            } else {
                sig.term <- dge.cols[dge.cols %in% c("padj", "FDR", "svalue", "adj.P.Val")]
                # If multiple matches, just use first
                if (length(sig.term) > 1) {
                    sig.term <- sig.term[1]
                }
            }
        }

        if (exp.term == "auto") {
            if (!any(dge.cols %in% c("baseMean", "logCPM", "AveExpr"))) {
                stop("Cannot determine significance term, please provide the column name to exp.term")
            } else {
                exp.term <- dge.cols[dge.cols %in% c("baseMean", "logCPM", "AveExpr")]
                # If multiple matches, just use first
                if (length(exp.term) > 1) {
                    exp.term <- exp.term[1]
                }
            }
        }
    }

    # Plot data.
    pd <- plotEnrichmentData(pathway = pathway.genes, stats = stats, gseaParam = gseaParam)

    rnk <- rank(-stats)
    ord <- order(rnk)

    statsAdj <- stats[ord]

    pathway.genes <- unname(as.vector(na.omit(match(pathway.genes, names(statsAdj)))))
    pathway.genes <- sort(pathway.genes)
    pathway.genes <- unique(pathway.genes)

    gene.ids <- names(statsAdj[pathway.genes])
    pd$ticks$gene <- gene.ids
    pd$ticks$color <- tick.color

    # Color by DE status if DE results are provided.
    if (!is.null(dge.res)) {
        if (id.term == "rownames") {
            dge.res$ID <- rownames(dge.res)
        } else {
            dge.res$ID <- dge.res[[id.term]]
        }

        dge.res <- dge.res[match(gene.ids, dge.res$ID), ]

        # Get up and downregulated genes.
        up <- dge.res$ID[dge.res[[lfc.term]] > lfc.thresh & dge.res[[sig.term]] < sig.thresh & dge.res[[exp.term]] > exp.thresh]
        down <- dge.res$ID[dge.res[[lfc.term]] < -lfc.thresh & dge.res[[sig.term]] < sig.thresh & dge.res[[exp.term]] > exp.thresh]

        # Color by DE status.
        pd$ticks$color <- ifelse(pd$ticks$gene %in% up, dge.up.color,
            ifelse(pd$ticks$gene %in% down, dge.down.color, tick.color)
        )
    }

    p <- with(
        pd,
        ggplot(data = curve) +
            geom_ribbon(
                data = stats,
                mapping = aes(
                    x = rank, ymin = 0,
                    ymax = stat / maxAbsStat * (spreadES / 4)
                ),
                fill = "grey"
            ) +
            geom_line(aes(x = rank, y = ES), color = "green") +
            geom_segment(
                data = ticks,
                mapping = aes(
                    x = rank, y = -spreadES / 16,
                    xend = rank, yend = spreadES / 16,
                    text = gene, color = color
                ),
                size = 0.2
            ) +
            scale_color_identity() +
            geom_hline(yintercept = posES, colour = "red", linetype = "dashed") +
            geom_hline(yintercept = negES, colour = "red", linetype = "dashed") +
            geom_hline(yintercept = 0, colour = "black") +
            theme(
                panel.background = element_blank(),
                panel.grid.major = element_line(color = "grey92")
            ) +
            labs(x = "Rank", y = "Enrichment Score", title = pathway.name)
    )

    # Add plot border, add ticks, set axis labels.
    ay <- list(
        showline = TRUE,
        mirror = TRUE,
        linecolor = toRGB("black"),
        linewidth = 0.5,
        showgrid = FALSE
    )

    ax <- list(
        showline = TRUE,
        mirror = TRUE,
        linecolor = toRGB("black"),
        linewidth = 0.5,
        showgrid = FALSE,
    )

    fig <- ggplotly(p, tooltip = c("x", "text")) %>%
        config(
            edits = list(
                annotationPosition = TRUE,
                annotationTail = TRUE
            ),
            toImageButtonOptions = list(format = "svg"),
            displaylogo = FALSE,
            plotGlPixelRatio = 7
        ) %>%
        layout(
            showlegend = FALSE,
            xaxis = ax,
            yaxis = ay
        )

    # Feature count annotations.
    if (!is.null(fgsea.res) && plot.stats) {
        padj <- fgsea.res$padj[fgsea.res$pathway == pathway.name]
        ES <- fgsea.res$ES[fgsea.res$pathway == pathway.name]
        NES <- fgsea.res$NES[fgsea.res$pathway == pathway.name]
        size <- fgsea.res$size[fgsea.res$pathway == pathway.name]

        if (ES < 0) {
            anno.x <- 0
            anno.y <- 0.05
        } else {
            anno.x <- 1
            anno.y <- 0.95
        }

        fig <- fig %>%
            add_annotations(
                x = anno.x,
                y = anno.y,
                xref = "paper",
                yref = "paper",
                text = paste0(
                    "padj: ", padj,
                    "\nES: ", ES,
                    "\nNES: ", NES,
                    "\nGeneset Size: ", size
                ),
                showarrow = FALSE,
                font = list(size = 10)
            )
    }

    fig
}


#' Create barplot for top significant terms from fgsea analysis
#'
#' This creates a useful, easy to interpret barplot to summarize the top
#' significant gene sets from an fgsea analysis.
#'
#' @param gsea.res A data.frame of GSEA results as returned by `runGSEA` or `runCustomGSEA`.
#' @param padj.th The significance threshold (adjusted p-value) for filtering gene sets.
#' @param top The number of top significant gene sets to consider.
#'
#' @return A plotly plot.
#'
#' @importFrom ggplot2 geom_col coord_flip labs theme_bw scale_fill_viridis ylim theme
#'   element_text scale_x_discrete
#' @importFrom plotly ggplotly
#' 
#' @author Jared Andrews
#' @export
plot_gsea_barplot <- function(gsea.res, padj.th = 0.05, top = 50) {
    df.sub <- gsea.res[gsea.res$padj < padj.th, ]

    if (nrow(df.sub) > top) {
        df.sub <- df.sub[order(df.sub$padj), ]
        df.sub <- df.sub[1:top, ]
    }

    if (nrow(df.sub) > 0) {

        p <- ggplot(df.sub, aes(reorder(pathway, -NES), NES)) +
            geom_col(aes(fill = -log10(padj), text = padj)) +
            coord_flip() +
            labs(
                x = NULL, y = "Normalized Enrichment Score",
                title = paste0(ct, " - Top ", top, "\np.adj < ", padj.th)
            ) +
            theme_bw() +
            scale_fill_viridis() +
            ylim(-5, 5) +
            theme(axis.text.y = element_text(size = 6), plot.title = element_text(size = 10)) +
            scale_x_discrete(label = function(x) strwrap(x, width = 40, exdent = 1))

        fig <- ggplotly(p, text = c("y", "x", "text"))
    } else {
        fig <- .empty_plot(paste0("No significant gene sets with padj < ", padj.th, "."))
    }

    fig
}
