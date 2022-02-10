# For adding various lines to plots easily.
.vline <- function(x = 0, color = "red", width = 1, dash = "solid") {
  list(
    type = "line",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    x0 = x,
    x1 = x,
    line = list(color = color, width = width, dash = dash)
  )
}

.hline <- function(y = 0, color = "blue", width = 1, dash = "solid") {
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color, width = width, dash = dash)
  )
}

.fitline <- function(df, color = "black", width = 0.75, dash = "solid") {
  list(
    type = "line",
    line = list(color = color, width = width, dash = dash),
    xref = "x",
    yref = "y",
    y0 = min(df$fv),
    y1 = max(df$fv),
    x0 = df$lfc.x[df$fv == min(df$fv)],
    x1 = df$lfc.x[df$fv == max(df$fv)]
  )
}

.make_volcano <- function(res, xlim, ylim, fc.thresh, fc.lines,
                          sig.line, h.id, feat.term, sig.term, lfc.term, down.color, up.color,
                          insig.color, sig.thresh = 0.05, fs = NULL, sig.size, insig.size,
                          sig.opacity, insig.opacity, label.size, webgl, webgl.ratio, show.counts,
                          show.hl.counts, counts.size, highlight.featsets, highlight.feats, featsets,
                          highlight.feats.color, highlight.feats.size, highlight.feats.opac,
                          highlight.feats.linecolor, highlight.feats.linewidth,
                          highlight.featsets.color, highlight.featsets.size, highlight.featsets.opac,
                          highlight.featsets.linecolor, highlight.featsets.linewidth) {

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
  highlight <- NULL
  if (!is.null(highlight.feats) & highlight.feats != "") {
    highlight.feats <- strsplit(highlight.feats, ",|\\s|,\\s")[[1]]
    highlight <- highlight.feats[highlight.feats != ""]
  }

  highlight.fs <- NULL
  if (!is.null(highlight.featsets)) {
    for (featset in highlight.featsets) {
      highlight.fs <- c(highlight.fs, featsets[[featset]])
    }
  }

  # Significance filter.
  up.feats <- res[[sig.term]] < sig.thresh & res[[lfc.term]] > 0
  res$col[up.feats] <- up.color
  res$cex[up.feats] <- sig.size
  res$order[up.feats] <- 1
  res$opacity[up.feats] <- sig.opacity

  dn.feats <- res[[sig.term]] < sig.thresh & res[[lfc.term]] < 0
  res$col[dn.feats] <- down.color
  res$cex[dn.feats] <- sig.size
  res$order[dn.feats] <- 1
  res$opacity[dn.feats] <- sig.opacity

  # LFC filter.
  if(fc.thresh > 0) {
    fc.threshed <- abs(res[[lfc.term]]) < fc.thresh
    res$col[fc.threshed] <- insig.color
    res$cex[fc.threshed] <- insig.size
    res$order[fc.threshed] <- 0
  }

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

  res$hover.string <- paste("</br><b>", feat.term, ":</b> ", res$feat,
                            "</br><b>", lfc.term, ":</b> ", format(round(res[[lfc.term]], 4), nsmall = 4),
                            "</br><b>", sig.term, ":</b> ", format(round(res[[sig.term]], 6), nsmall = 6))
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
    sig.hline <- .hline(y = -log10(sig.thresh), color = "#999999", width = 1, dash = "longdash")
  }

  if (fc.thresh != 0 & fc.lines) {
    fc.line1 <- .vline(x = fc.thresh, color = "#999999", width = 1, dash = "longdash")
    fc.line2 <- .vline(x = -fc.thresh, color = "#999999", width = 1, dash = "longdash")
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
                 source = paste0(h.id, "_volc")) %>%
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
             shapes = list(sig.hline, fc.line1, fc.line2)) %>%
      add_annotations(x = fs$x, y = fs$y, text = gs$customdata,
                      font = list(size = label.size, family = "Arial"), arrowside = "none")
  } else {
    fig <- fig %>% layout(xaxis = ax,
                          yaxis = ay,
                          showlegend = FALSE,
                          shapes = list(sig.hline, fc.line1, fc.line2))
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
