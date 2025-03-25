# Make the heatmap for differentially expressed genes under certain cutoffs.
.make_heatmap <- function(mat, res, anno, bm.col.func, lfc.col.func, sig.term,
                          sig.thresh = 0.05, base_mean = 0, log2fc = 0, row.km = 0, col.km = 0) {

  # Adjust for potential differences in the results table.
  l <- res[[sig.term]] <= sig.thresh & res$baseMean >= base_mean & abs(res$log2FoldChange) >= log2fc; l[is.na(l)] = FALSE


  # If 0 or 1 gene meet the cutoffs, don't make heatmap.
  if(sum(l) < 2) {
    return(NULL)
  }

  m <- mat[l, ]
  m <- t(scale(t(m)))

  env$row_index <- which(l)

  # Sets color palette to dittoSeq colors instead of random
  if (!is.null(anno)) {
    ds.colors <- dittoColors()
    anno.colors <- list()
    i <- 1

    for (n in names(anno)) {
      out <- list()

      for (lev in unique(anno[[n]])) {
        out[[lev]] <- ds.colors[i]
        i <- i + 1
      }

      anno.colors[[n]] <- unlist(out)
    }

    anno <- HeatmapAnnotation(df = anno, col = anno.colors)
  }

  basem_df <- log10(res$baseMean[l]+1)
  names(basem_df) <- rownames(m)

  lfc_df <- res$log2FoldChange[l]
  names(lfc_df) <- rownames(m)

  ht <- Heatmap(m, name = "z-score",
                top_annotation = anno,
                show_row_names = FALSE, show_column_names = FALSE,
                row_km = row.km, column_km = col.km,
                column_title_gp = gpar(fontsize = 10),
                column_title = paste0(sum(l), " significant genes \nwith ", sig.term," < ", sig.thresh),
                show_row_dend = FALSE) +
    Heatmap(basem_df, show_row_names = FALSE, width = unit(5, "mm"),
            name = "log10(baseMean+1)", col = bm.col.func, show_column_names = FALSE) +
    Heatmap(lfc_df, show_row_names = FALSE, width = unit(5, "mm"),
            name = "log2FoldChange", col = lfc.col.func, show_column_names = FALSE)
  ht <- draw(ht, merge_legend = TRUE)
  ht
}


# Generic volcano plot function.
# TODO: Add defaults, document, and export.
.make_volcano <- function(res, xlim, ylim, fc.thresh, fc.lines, hover.info = NULL,
                          sig.line, h.id, feat.term, sig.term, lfc.term, down.color, up.color,
                          insig.color, sig.thresh = 0.05, basemean.thresh = 0, fs = NULL, sig.size, insig.size,
                          sig.opacity, insig.opacity, label.size, webgl, webgl.ratio, show.counts,
                          show.hl.counts, counts.size, highlight.featsets, highlight.feats, featsets,
                          highlight.feats.color, highlight.feats.size, highlight.feats.opac,
                          highlight.feats.linecolor, highlight.feats.linewidth,
                          highlight.featsets.color, highlight.featsets.size, highlight.featsets.opac,
                          highlight.featsets.linecolor, highlight.featsets.linewidth, h.id.suffix = "_volc") {

  # Styling.
  res$col <- rep(insig.color, nrow(res))
  res$cex <- rep(insig.size, nrow(res))
  res$order <- rep(0, nrow(res))
  res$lcol <- res$col
  res$lw <- 0
  res$opacity <- insig.opacity

  # Remove features with NA significance term (due to low expression, etc).
  res <- res[!is.na(res[[sig.term]]),]

  # Get all gene IDs or symbols to be highlighted.
  highlight <- highlight.feats
  if (!is.null(highlight.feats) & highlight.feats != "") {
    highlight.feats <- strsplit(highlight.feats, ",|\\s|,\\s")[[1]]
    highlight <- highlight.feats[highlight.feats != ""]
  }

  highlight.fs <- highlight.featsets
  if (!is.null(highlight.featsets)) {
    for (featset in highlight.featsets) {
      highlight.fs <- c(highlight.fs, featsets[[featset]])
    }
  }

  # Significance filter.
  up.feats <- res[[sig.term]] < sig.thresh & res[[lfc.term]] > fc.thresh & res$baseMean > basemean.thresh
  res$col[up.feats] <- up.color
  res$cex[up.feats] <- sig.size
  res$order[up.feats] <- 1
  res$opacity[up.feats] <- sig.opacity

  dn.feats <- res[[sig.term]] < sig.thresh & res[[lfc.term]] < -fc.thresh & res$baseMean > basemean.thresh
  res$col[dn.feats] <- down.color
  res$cex[dn.feats] <- sig.size
  res$order[dn.feats] <- 1
  res$opacity[dn.feats] <- sig.opacity

  res$x <- res[[lfc.term]]
  res$y <- -log10(res[[sig.term]])

  res$col[res$y < -log10(sig.thresh)] <- insig.color

  res$sh <- ifelse(res$y > ylim, "triangle-up-open",
                   ifelse(res$x < -xlim, "triangle-left-open",
                          ifelse(res$x > xlim, "triangle-right-open", 0)))

  res$lw <- ifelse(res$sh != 0, 1, 0)

  res$y[res$y > ylim] <- ylim - 0.2
  res$x[res$x > xlim] <- xlim - 0.05
  res$x[res$x < -xlim] <- -xlim + 0.05
  if (feat.term == "rows") {
    res$feat <- rownames(res)
  } else {
    res$feat <- res[[feat.term]]
  }

  # Gene/geneset highlighting.
  n.fs.hl <- 0
  n.hl <- 0

  if (!is.null(highlight.fs)) {
    highlight.fs <- highlight.fs[highlight.fs %in% res$feat]
    n.fs.hl <- length(res$col[res$feat %in% highlight.fs])

    res$col[res$feat %in% highlight.fs] <- highlight.featsets.color
    res$cex[res$feat %in% highlight.fs] <- highlight.featsets.size
    res$opacity[res$feat %in% highlight.fs] <- highlight.featsets.opac
    res$lcol[res$feat %in% highlight.fs] <- highlight.featsets.linecolor
    res$lw[res$feat %in% highlight.fs] <- highlight.featsets.linewidth
    res$order[res$feat %in% highlight.fs] <- 2
  }

  # Want these to have precedence over the feature sets in case entries are in both.
  if (!is.null(highlight)) {
    highlight <- highlight[highlight %in% res$feat]
    n.hl <- length(res$col[res$feat %in% highlight])

    res$col[res$feat %in% highlight] <- highlight.feats.color
    res$cex[res$feat %in% highlight] <- highlight.feats.size
    res$opacity[res$feat %in% highlight] <- highlight.feats.opac
    res$lcol[res$feat %in% highlight] <- highlight.feats.linecolor
    res$lw[res$feat %in% highlight] <- highlight.feats.linewidth
    res$order[res$feat %in% highlight] <- 3
  }

  res$hover.string <- paste0("</br><b>", feat.term, ":</b> ", res$feat,
                             "</br><b>", lfc.term, ":</b> ", format(round(res[[lfc.term]], 4), nsmall = 4),
                             "</br><b>", sig.term, ":</b> ", format(round(res[[sig.term]], 6), nsmall = 6))

  if (!is.null(hover.info)) {
    for (n in hover.info) {
      res$hover.string <- paste0(res$hover.string, "</br><b>", n, ":</b>", res[[n]])
    }
  }

  res <- as.data.frame(res)
  res <- res[order(res$order),]

  # Get feature numbers.
  n.up.feats <- length(up.feats[up.feats == TRUE])
  n.dn.feats <- length(dn.feats[dn.feats == TRUE])
  n.feats <- nrow(res)

  # Add plot border.
  ay <- list(
    showline = TRUE,
    mirror = TRUE,
    linecolor = toRGB("black"),
    linewidth = 0.5,
    title = paste0("-log10(", sig.term, ")"),
    range = list(0, ylim),
    showgrid = FALSE,
    layer = "below traces",
    zeroline = FALSE,
    ticks = "outside",
    zerolinewidth = 0.5
  )

  ax <- list(
    showline = TRUE,
    mirror = TRUE,
    linecolor = toRGB("black"),
    linewidth = 0.5,
    title = lfc.term,
    range = list(-xlim, xlim),
    showgrid = FALSE,
    layer = "below traces",
    ticks = "outside",
    zerolinewidth = 0.5
  )

  # Create vertical and horizontal lines.
  fc.line1 <- NULL
  fc.line2 <- NULL

  sig.hline <- NULL
  if(sig.line) {
    sig.hline <- .hline(y = -log10(sig.thresh), color = "#000000", width = 1, dash = "longdash")
  }

  if (fc.thresh != 0 & fc.lines) {
    fc.line1 <- .vline(x = fc.thresh, color = "#000000", width = 1, dash = "longdash")
    fc.line2 <- .vline(x = -fc.thresh, color = "#000000", width = 1, dash = "longdash")
  }

  # Figure generation.
  fig <- plot_ly(res, x = ~x,
                 y = ~y,
                 customdata = ~feat,
                 type = "scatter",
                 mode = "markers",
                 marker = list(color = ~col,
                               size = ~cex,
                               symbol = ~sh,
                               line = list(color = ~lcol, width = ~lw),
                               opacity = ~opacity),
                 text = ~hover.string,
                 hoverinfo = "text",
                 source = paste0(h.id, h.id.suffix)) %>%
    config(edits = list(annotationPosition = TRUE,
                        annotationTail = TRUE),
           toImageButtonOptions = list(format = "svg"),
           displaylogo = FALSE,
           plotGlPixelRatio = webgl.ratio)

  if (!is.null(fs)) {
    fig <- fig %>%
      layout(xaxis = ax,
             yaxis = ay,
             showlegend = FALSE,
             shapes = list(sig.hline, fc.line1, fc.line2),
             hoverlabel = list(font=list(size=10))) %>%
      add_annotations(x = fs$x, y = fs$y, text = fs$customdata,
                      font = list(size = label.size, family = "Arial"), arrowside = "none")
  } else {
    fig <- fig %>% layout(xaxis = ax,
                          yaxis = ay,
                          showlegend = FALSE,
                          shapes = list(sig.hline, fc.line1, fc.line2),
                          hoverlabel = list(font=list(size=10)))
  }

  # Feature count annotations.
  if (show.counts) {
    fig <- fig %>%
      add_annotations(
        x= 1,
        y= 1,
        xref = "paper",
        yref = "paper",
        text = paste0("Up features: ", n.up.feats,
                      "\nDown features: ", n.dn.feats,
                      "\nTotal features: ", n.feats),
        showarrow = FALSE,
        font = list(size = counts.size)
      )
  }

  if (show.hl.counts) {
    fig <- fig %>%
      add_annotations(
        x= 0,
        y= 1,
        xref = "paper",
        yref = "paper",
        text = paste0("Set features: ", n.fs.hl,
                      "\nHighlighted features: ", n.hl),
        showarrow = FALSE,
        font = list(size = counts.size)
      )
  }

  if (webgl) {
    fig <- fig %>% toWebGL()
  }

  fig
}


.make_maplot <- function(res, ylim, fc.thresh, fc.lines, h.id, sig.term,
                         down.color, up.color, insig.color, sig.size, insig.size, sig.thresh = 0.05,
                         gs = NULL, basemean.thresh = 0, sig.opacity, insig.opacity, label.size, webgl, webgl.ratio, show.counts,
                         counts.size, show.hl.counts, highlight.genesets, highlight.genes, genesets,
                         highlight.genes.color, highlight.genes.size, highlight.genes.opac,
                         highlight.genes.linecolor, highlight.genes.linewidth,
                         highlight.genesets.color, highlight.genesets.size, highlight.genesets.opac,
                         highlight.genesets.linecolor, highlight.genesets.linewidth,
                         loess, loess.color, loess.span, loess.genesets, loess.genesets.color, loess.genesets.span) {

  res$col <- rep(insig.color, nrow(res))
  res$cex <- rep(insig.size, nrow(res))
  res$order <- rep(0, nrow(res))
  res$lcol <- res$col
  res$lw <- 0
  res$opacity <- insig.opacity

  # Remove genes with NA padj/svalue due to low expression.
  res <- res[!is.na(res[[sig.term]]),]

  # Get all gene IDs or symbols to be highlighted.
  highlight <- NULL
  if (!is.null(highlight.genes) & highlight.genes != "") {
    highlight.genes <- strsplit(highlight.genes, ",|\\s|,\\s")[[1]]
    highlight <- highlight.genes[highlight.genes != ""]
  }

  highlight.gs <- NULL
  if (!is.null(highlight.genesets)) {
    for (geneset in highlight.genesets) {
      highlight.gs <- c(highlight.gs, genesets[[geneset]])
    }
  }

  # Styling.
  up.degs <- res[[sig.term]] < sig.thresh & res$log2FoldChange > fc.thresh & res$baseMean > basemean.thresh
  res$col[up.degs] <- up.color
  res$cex[up.degs] <- sig.size
  res$order[up.degs] <- 1
  res$opacity[up.degs] <- sig.opacity

  dn.degs <- res[[sig.term]] < sig.thresh & res$log2FoldChange < -fc.thresh & res$baseMean > basemean.thresh
  res$col[dn.degs] <- down.color
  res$cex[dn.degs] <- sig.size
  res$order[dn.degs] <- 1
  res$opacity[dn.degs] <- sig.opacity

  res$lcol <- res$col

  # Get gene numbers.
  n.up.genes <- length(up.degs[up.degs == TRUE])
  n.dn.genes <- length(dn.degs[dn.degs == TRUE])
  n.genes <- nrow(res)

  res$x <- res$baseMean
  res$y <- res$log2FoldChange
  res$sh <- ifelse(res$log2FoldChange > ylim, "triangle-up-open",
                   ifelse(res$log2FoldChange < -ylim, "triangle-down-open", 0))
  res$lw <- ifelse(res$sh != 0, 1, 0)
  res$y[res$y > ylim] <- ylim - 0.05
  res$y[res$y < -ylim] <- -ylim + 0.05
  res$Gene <- rownames(res)

  # Gene/geneset highlighting.
  n.gs.hl <- 0
  n.hl <- 0

  highlight.data <- NULL
  
  if (!is.null(highlight.gs)) {
    highlight.gs <- highlight.gs[highlight.gs %in% res$Gene]
    n.gs.hl <- length(res$col[res$Gene %in% highlight.gs])

    res$col[res$Gene %in% highlight.gs] <- highlight.genesets.color
    res$cex[res$Gene %in% highlight.gs] <- highlight.genesets.size
    res$opacity[res$Gene %in% highlight.gs] <- highlight.genesets.opac
    res$lcol[res$Gene %in% highlight.gs] <- highlight.genesets.linecolor
    res$lw[res$Gene %in% highlight.gs] <- highlight.genesets.linewidth
    res$order[res$Gene %in% highlight.gs] <- 2

    highlight.data <- res[res$Gene %in% highlight.gs,]
  }

  # Want these to have precedence over the genesets in case entries are in both.
  if (!is.null(highlight)) {
    highlight <- highlight[highlight %in% res$Gene]
    n.hl <- length(res$col[res$Gene %in% highlight])

    res$col[res$Gene %in% highlight] <- highlight.genes.color
    res$cex[res$Gene %in% highlight] <- highlight.genes.size
    res$opacity[res$Gene %in% highlight] <- highlight.genes.opac
    res$lcol[res$Gene %in% highlight] <- highlight.genes.linecolor
    res$lw[res$Gene %in% highlight] <- highlight.genes.linewidth
    res$order[res$Gene %in% highlight] <- 3
  }

  res$hover.string <- paste("</br><b>Gene:</b> ", res$Gene,
                            "</br><b>log2 Fold Change:</b> ", format(round(res$log2FoldChange, 4), nsmall = 4),
                            "</br><b>", sig.term, ":</b> ", format(round(res[[sig.term]], 4), nsmall = 4),
                            "</br><b>baseMean (avg. Expression):</b> ", format(round(res$baseMean, 2), nsmall = 2))

  res <- as.data.frame(res)
  res <- res[order(res$order),]

  # Add plot border.
  ay <- list(
    showline = TRUE,
    mirror = TRUE,
    linecolor = toRGB("black"),
    linewidth = 0.5,
    title = "log2(fold change)",
    range = list(-ylim, ylim),
    showgrid = FALSE,
    layer = "below traces",
    ticks = "outside",
    zerolinewidth = 0.5
  )

  ax <- list(
    showline = TRUE,
    mirror = TRUE,
    linecolor = toRGB("black"),
    linewidth = 0.5,
    title = "log10(baseMean)",
    showgrid = FALSE,
    layer = "below traces",
    zeroline = FALSE,
    ticks = "outside",
    zerolinewidth = 0.5
  )

  # Create vertical and horizontal lines.
  fc.line1 <- NULL
  fc.line2 <- NULL
  if (fc.thresh != 0 & fc.lines) {
    fc.line1 <- .hline(y = fc.thresh, color = "#999999", width = 1, dash = "longdash")
    fc.line2 <- .hline(y = -fc.thresh, color = "#999999", width = 1, dash = "longdash")
  }
  
  fig <- plot_ly(res, x = ~log10(x),
                 source = paste0(h.id, "_ma")) %>%
                 add_markers(y = ~y,
                 customdata = ~Gene,
                 marker = list(color = ~col,
                               size = ~cex,
                               symbol = ~sh,
                               line = list(color = ~lcol, width = ~lw),
                               opacity = ~opacity),
                 text = ~hover.string,
                 hoverinfo = "text")


  if (loess) {
    fig <- fig %>% add_lines(y = ~fitted(loess(y ~ -log10(x), span = loess.span)), line = list(color = loess.color, width = 3), name = "LOESS")
  }

  if (!is.null(highlight.data)) {
    if (loess.genesets) {
      fig <- fig %>% add_lines(data = highlight.data, x = ~log10(x), y = ~fitted(loess(y ~ -log10(x), span = loess.genesets.span)), 
        line = list(color = loess.genesets.color, width = 3), name = "LOESS")
    }
  }

  fig <- fig %>% 
    config(edits = list(annotationPosition = TRUE,
                        annotationTail = TRUE),
           toImageButtonOptions = list(format = "svg"),
           displaylogo = FALSE,
           plotGlPixelRatio = webgl.ratio)

  if (!is.null(gs)) {
    fig <- fig %>%
      layout(xaxis = ax,
             yaxis = ay,
             showlegend = FALSE,
             shapes = list(fc.line1, fc.line2),
             hoverlabel = list(font=list(size=10))) %>%
      add_annotations(x = gs$x, y = gs$y, text = gs$customdata,
                      font = list(size = label.size, family = "Arial"), arrowside = "none")
  } else {
    fig <- fig %>% layout(xaxis = ax,
                   yaxis = ay,
                   showlegend = FALSE,
                   shapes = list(fc.line1, fc.line2),
                   hoverlabel = list(font=list(size=10)))
  }

  # Gene count annotations.
  if (show.counts) {
    fig <- fig %>%
      add_annotations(
        x= 1,
        y= 0,
        xref = "paper",
        yref = "paper",
        text = paste0("Up-reg. genes: ", n.up.genes,
                      "\nDown-reg. genes: ", n.dn.genes,
                      "\nTotal genes: ", n.genes),
        showarrow = FALSE,
        font = list(size = counts.size)
      )
  }

  if (show.hl.counts) {
    fig <- fig %>%
      add_annotations(
        x= 1,
        y= 1,
        xref = "paper",
        yref = "paper",
        text = paste0("Geneset genes: ", n.gs.hl,
                      "\nHighlighted genes: ", n.hl),
        showarrow = FALSE,
        font = list(size = counts.size)
      )
  }

  if (webgl) {
    fig <- fig %>% toWebGL()
  }

  fig
}


# Determine columns in results table and adjust output appropriately.
.make_gene_output <- function(df) {
  if("svalue" %in% colnames(df)) {
    tmpl <- "
<pre>
Gene: @{gene}
baseMean: @{df[1, 'baseMean']}
log2FoldChange: @{df[1, 'log2FoldChange']}
lfcSE: @{df[1, 'lfcSE']}
svalue: @{df[1, 'svalue']}</pre>
"
  } else if(!("stat" %in% colnames(df))) {
    tmpl <- "
<pre>
Gene: @{gene}
baseMean: @{df[1, 'baseMean']}
log2FoldChange: @{df[1, 'log2FoldChange']}
lfcSE: @{df[1, 'lfcSE']}
pvalue: @{df[1, 'pvalue']}
adjusted pvalue: @{df[1, 'padj']}</pre>
"
  } else {
    tmpl <- "
<pre>
Gene: @{gene}
baseMean: @{df[1, 'baseMean']}
log2FoldChange: @{df[1, 'log2FoldChange']}
lfcSE: @{df[1, 'lfcSE']}
stat: @{df[1, 'stat']}
pvalue: @{df[1, 'pvalue']}
adjusted pvalue: @{df[1, 'padj']}</pre>
"
  }

  return(tmpl)
}
