# Make the heatmap for differentially expressed genes under certain cutoffs.
.make_heatmap <- function(mat, res, anno, bm.col.func, lfc.col.func,
                          fdr = 0.05, base_mean = 0, log2fc = 0, row.km = 0, col.km = 0) {

  # Adjust for potential differences in the results table.
  sig.term <- "padj"
  if("svalue" %in% colnames(res)) {
    l <- res$svalue <= fdr & res$baseMean >= base_mean & abs(res$log2FoldChange) >= log2fc; l[is.na(l)] = FALSE
    sig.term <- "svalue"
  } else {
    l <- res$padj <= fdr & res$baseMean >= base_mean & abs(res$log2FoldChange) >= log2fc; l[is.na(l)] = FALSE
  }

  # If no genes meet the cutoffs, don't make heatmap.
  if(sum(l) == 0) {
    return(NULL)
  }

  m <- mat[l, ]
  m <- t(scale(t(m)))

  env$row_index <- which(l)

  ht <- Heatmap(m, name = "z-score",
                top_annotation = HeatmapAnnotation(df = anno),
                show_row_names = FALSE, show_column_names = FALSE,
                row_km = row.km, column_km = col.km,
                column_title_gp = gpar(fontsize = 10),
                column_title = paste0(sum(l), " significant genes \nwith ", sig.term," < ", fdr),
                show_row_dend = FALSE) +
    Heatmap(log10(res$baseMean[l]+1), show_row_names = FALSE, width = unit(5, "mm"),
            name = "log10(baseMean+1)", col = bm.col.func, show_column_names = FALSE) +
    Heatmap(res$log2FoldChange[l], show_row_names = FALSE, width = unit(5, "mm"),
            name = "log2FoldChange", col = lfc.col.func, show_column_names = FALSE)
  ht <- draw(ht, merge_legend = TRUE)
  ht
}


.make_maplot <- function(res, ylim, fc.thresh, fc.lines, h.id, sig.term,
                         down.color, up.color, insig.color, sig.thresh = 0.05, gs = NULL) {

  res$col <- rep(insig.color, nrow(res))
  res$cex <- rep(3, nrow(res))
  res$order <- rep(0, nrow(res))

  # Significance filter.
  up.degs <- res[[sig.term]] < sig.thresh & res$log2FoldChange > 0
  res$col[up.degs] <- up.color
  res$cex[up.degs] <- 5
  res$order[up.degs] <- 1

  dn.degs <- res[[sig.term]] < sig.thresh & res$log2FoldChange < 0
  res$col[dn.degs] <- down.color
  res$cex[dn.degs] <- 5
  res$order[dn.degs] <- 1

  # LFC filter.
  if(fc.thresh > 0) {
    fc.threshed <- abs(res$log2FoldChange) < fc.thresh
    res$col[fc.threshed] <- insig.color
    res$cex[fc.threshed] <- 3
    res$order[fc.threshed] <- 0
  }

  # Remove genes with NA padj/svalue due to low expression.
  res <- res[!is.na(res[[sig.term]]),]

  res$x <- res$baseMean
  res$y <- res$log2FoldChange
  res$sh <- ifelse(res$log2FoldChange > ylim, "triangle-up-open",
                   ifelse(res$log2FoldChange < -ylim, "triangle-down-open", 0))
  res$y[res$y > ylim] <- ylim - 0.05
  res$y[res$y < -ylim] <- -ylim + 0.05
  res$Gene <- rownames(res)

  res$hover.string <- paste("</br><b>Gene:</b> ", res$Gene,
                            "</br><b>log2 Fold Change:</b> ", format(round(res$log2FoldChange, 4), nsmall = 4),
                            "</br><b>", sig.term, ":</b> ", format(round(res[[sig.term]], 4), nsmall = 4),
                            "</br><b>baseMean (avg. Expression):</b> ", format(round(res$baseMean, 2), nsmall = 2))

  res <- as.data.frame(res)
  res <- res[order(res$order),]

  # Add plot border.
  ay <- list(
    showline = TRUE,
    mirror = "ticks",
    linecolor = toRGB("black"),
    linewidth = 1,
    title = "log2(fold change)",
    range = list(-ylim, ylim),
    showgrid = FALSE,
    layer = "below traces"
  )

  ax <- list(
    showline = TRUE,
    mirror = "ticks",
    linecolor = toRGB("black"),
    linewidth = 1,
    title = "log10(baseMean)",
    showgrid = FALSE,
    layer = "below traces",
    zeroline = FALSE
  )

  # Create vertical and horizontal lines.
  fc.line1 <- NULL
  fc.line2 <- NULL
  if (fc.thresh != 0 & fc.lines) {
    fc.line1 <- .hline(y = fc.thresh, color = "#999999", width = 1, dash = "longdash")
    fc.line2 <- .hline(y = -fc.thresh, color = "#999999", width = 1, dash = "longdash")
  }

  fig <- plot_ly(res, x = ~log10(x),
                 y = ~y,
                 customdata = ~Gene,
                 type = "scatter",
                 mode = "markers",
                 marker = list(color = ~col,
                               size = ~cex,
                               symbol = ~sh,
                               line = list(color = ~col)),
                 text = ~hover.string,
                 hoverinfo = "text",
                 source = paste0(h.id, "_ma")) %>%
    config(edits = list(annotationPosition = TRUE,
                        annotationTail = TRUE),
           toImageButtonOptions = list(format = "svg"))

  if (!is.null(gs)) {
    fig %>%
      layout(xaxis = ax,
             yaxis = ay,
             showlegend = FALSE,
             shapes = list(fc.line1, fc.line2)) %>%
      add_annotations(x = gs$x, y = gs$y, text = gs$customdata,
                      font = list(size = 10, family = "Arial"), arrowside = "none") %>%
      toWebGL()
  } else {
    fig %>% layout(xaxis = ax,
                   yaxis = ay,
                   showlegend = FALSE,
                   shapes = list(fc.line1, fc.line2)) %>% toWebGL()
  }
}


# make the volcano plot with some genes highlighted
.make_volcano <- function(res, xlim, ylim, fc.thresh, fc.lines,
                          sig.line, h.id, sig.term, down.color, up.color,
                          insig.color, sig.thresh = 0.05, gs = NULL) {

  # Adjust for potential differences in the results table.
  res$col <- rep(insig.color, nrow(res))
  res$cex <- rep(3, nrow(res))
  res$order <- rep(0, nrow(res))

  # Significance filter.
  up.degs <- res[[sig.term]] < sig.thresh & res$log2FoldChange > 0
  res$col[up.degs] <- up.color
  res$cex[up.degs] <- 5
  res$order[up.degs] <- 1

  dn.degs <- res[[sig.term]] < sig.thresh & res$log2FoldChange < 0
  res$col[dn.degs] <- down.color
  res$cex[dn.degs] <- 5
  res$order[dn.degs] <- 1

  # LFC filter.
  if(fc.thresh > 0) {
    fc.threshed <- abs(res$log2FoldChange) < fc.thresh
    res$col[fc.threshed] <- insig.color
    res$cex[fc.threshed] <- 3
    res$order[fc.threshed] <- 0
  }

  # Remove genes with NA padj/svalue due to low expression.
  res <- res[!is.na(res[[sig.term]]),]

  res$x <- res$log2FoldChange
  res$y <- -log10(res[[sig.term]])

  res$col[res$y < -log10(sig.thresh)] <- insig.color

  res$sh <- ifelse(res$y > ylim, "triangle-up-open",
               ifelse(res$x < -xlim, "triangle-left-open",
                      ifelse(res$x > xlim, "triangle-right-open", 0)))

  res$y[res$y > ylim] <- ylim - 0.2
  res$x[res$x > xlim] <- xlim - 0.05
  res$x[res$x < -xlim] <- -xlim + 0.05
  res$Gene <- rownames(res)

  res$hover.string <- paste("</br><b>Gene:</b> ", res$Gene,
                            "</br><b>log2 Fold Change:</b> ", format(round(res$log2FoldChange, 4), nsmall = 4),
                            "</br><b>", sig.term, ":</b> ", format(round(res[[sig.term]], 6), nsmall = 6),
                            "</br><b>baseMean (avg. Expression):</b> ", format(round(res$baseMean, 2), nsmall = 2))
  res <- as.data.frame(res)
  res <- res[order(res$order),]

  # Add plot border.
  ay <- list(
    showline = TRUE,
    mirror = "ticks",
    linecolor = toRGB("black"),
    linewidth = 1,
    title = paste0("-log10(", sig.term, ")"),
    range = list(0, ylim),
    showgrid = FALSE,
    layer = "below traces",
    zeroline = FALSE
  )

  ax <- list(
    showline = TRUE,
    mirror = "ticks",
    linecolor = toRGB("black"),
    linewidth = 1,
    title = "log2(fold change)",
    range = list(-xlim, xlim),
    showgrid = FALSE,
    layer = "below traces"
  )

  # Create vertical and horizontal lines.
  fc.line1 <- NULL
  fc.line2 <- NULL

  sig.hline <- NULL
  if(sig.line) {
    sig.hline <- .hline(y = -log10(sig.thresh), color = "#999999", width = 1, dash = "longdash")
  }

  if (fc.thresh != 0 & fc.lines) {
    fc.line1 <- .vline(x = fc.thresh, color = "#999999", width = 1, dash = "longdash")
    fc.line2 <- .vline(x = -fc.thresh, color = "#999999", width = 1, dash = "longdash")
  }

  fig <- plot_ly(res, x = ~x,
                 y = ~y,
                 customdata = ~Gene,
                 type = "scatter",
                 mode = "markers",
                 marker = list(color = ~col,
                               size = ~cex,
                               symbol = ~sh,
                               line = list(color = ~col)),
                 text = ~hover.string,
                 hoverinfo = "text",
                 source = paste0(h.id, "_volc")) %>%
    config(edits = list(annotationPosition = TRUE,
                        annotationTail = TRUE),
           toImageButtonOptions = list(format = "svg"))

  if (!is.null(gs)) {
    fig %>%
      layout(xaxis = ax,
             yaxis = ay,
             showlegend = FALSE,
             shapes = list(sig.hline, fc.line1, fc.line2)) %>%
      add_annotations(x = gs$x, y = gs$y, text = gs$customdata,
                      font = list(size = 10, family = "Arial"), arrowside = "none") %>%
      toWebGL()
  } else {
    fig %>% layout(xaxis = ax,
                   yaxis = ay,
                   showlegend = FALSE,
                   shapes = list(sig.hline, fc.line1, fc.line2)) %>% toWebGL()
  }
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
