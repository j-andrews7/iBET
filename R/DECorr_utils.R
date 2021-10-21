.make_xyplot <- function(res1, res2, res1.color, res2.color, both.color, insig.color,
                         sig.thresh, lfc.thresh, gene.col, opacity, df.vars,
                         regr = TRUE, genes.labeled = NULL, ylim, xlim, show,
                         source, label.size, webgl) {

  comp1.name <- names(res1)
  comp1 <- res1[[1]]

  comp2.name <- names(res2)
  comp2 <- res2[[1]]

  # If gene column not defined, set to rownames.
  if (is.null(gene.col)) {
    comp1$Gene <- rownames(comp1)
    comp2$Gene <- rownames(comp2)
  } else {
    comp1$Gene <- comp1[[gene.col]]
    comp2$Gene <- comp2[[gene.col]]
  }

  # Make column grabbing easier later.
  comp1$lfc.x <- comp1[[df.vars$res1.lfc.col]]
  comp2$lfc.y <- comp2[[df.vars$res2.lfc.col]]

  comp1$sig.x <- comp1[[df.vars$res1.sig.col]]
  comp2$sig.y <- comp2[[df.vars$res2.sig.col]]

  comp1$exp.x <- comp1[[df.vars$res1.expr.col]]
  comp2$exp.y <- comp2[[df.vars$res2.expr.col]]

  # Remove NAs.
  comp1 <- comp1[!is.na(comp1$lfc.x) & !is.na(comp1$sig.x),]
  comp2 <- comp2[!is.na(comp2$lfc.y) & !is.na(comp2$sig.y),]

  # Create final plotting df.
  full.df <- merge(comp1, comp2, by = "Gene")

  # Set significance status.
  full.df$Sig <- "Not Significant"
  if (is.null(lfc.thresh)) {
    lfc.thresh <- 0
  }

  full.df$Sig <- ifelse(full.df$sig.x < sig.thresh & abs(full.df$lfc.x) >= lfc.thresh,
                        ifelse(full.df$sig.y < sig.thresh & abs(full.df$lfc.y) >= lfc.thresh,
                               "Both Significant", paste0(comp1.name, " Significant")),
                        ifelse(full.df$sig.y < sig.thresh & abs(full.df$lfc.y) >= lfc.thresh,
                               paste0(comp2.name, " Significant"), "Not Significant"))

  # Drop not significant genes if needed.
  if (!("Both Significant" %in% show)) {
    full.df <- full.df[full.df$Sig != "Both Significant",]
  }

  if (!("Not Significant" %in% show)) {
    full.df <- full.df[full.df$Sig != "Not Significant",]
  }

  if (!("X-axis Significant" %in% show)) {
    full.df <- full.df[full.df$Sig != paste0(comp1.name, " Significant"),]
  }

  if (!("Y-axis Significant" %in% show)) {
    full.df <- full.df[full.df$Sig != paste0(comp2.name, " Significant"),]
  }

  full.df$col <- rep(insig.color, nrow(full.df))
  full.df$cex <- rep(3, nrow(full.df))
  full.df$order <- rep(0, nrow(full.df))

  # Styling.
  full.df$col[full.df$Sig == paste0(comp1.name, " Significant")] <- res1.color
  full.df$col[full.df$Sig == paste0(comp2.name, " Significant")] <- res2.color
  full.df$col[full.df$Sig == "Both Significant"] <- both.color
  full.df$cex[full.df$Sig != "Not Significant"] <- 5
  full.df$order[full.df$Sig != "Not Significant"] <- 1
  full.df$order[full.df$Sig == "Both Significant"] <- 2

  full.df$sh <- ifelse(full.df$lfc.y > ylim, "triangle-up-open",
                       ifelse(full.df$lfc.y < -ylim, "triangle-down-open",
                          ifelse(full.df$lfc.x < -xlim, "triangle-left-open",
                            ifelse(full.df$lfc.x > xlim, "triangle-right-open", 0))))

  # Calculate regression line if needed, prior to change values due to axis limits.
  regr.line <- NULL
  if (regr) {
    full.df$fv <- lm(lfc.y ~ lfc.x, data = full.df) %>% fitted.values()
    regr.line <- .fitline(full.df, width = 0.5)

    regr.anno <- paste0("R = ",
                  round(with(full.df, cor.test(lfc.x, lfc.y))$estimate, 2),
                  "p = ",
                  format(with(full.df, cor.test(lfc.x, lfc.y))$p.value, scientific = TRUE, digits = 3))
  }

  full.df$lfc.y[full.df$lfc.y > ylim] <- ylim - 0.1
  full.df$lfc.y[full.df$lfc.y < -ylim] <- -ylim + 0.1

  full.df$lfc.x[full.df$lfc.x > xlim] <- xlim - 0.1
  full.df$lfc.x[full.df$lfc.x < -xlim] <- -xlim + 0.1

  full.df$hover.string <- paste("</br><b>Gene:</b> ", full.df$Gene,
                                "</br><b>", paste0(comp1.name, " ", df.vars$res1.lfc.col), ":</b> ",
                                format(round(full.df$lfc.x, 4), nsmall = 4, scientific = FALSE),
                                "</br><b>", paste0(comp2.name, " ", df.vars$res2.lfc.col), ":</b> ",
                                format(round(full.df$lfc.y, 4), nsmall = 4, scientific = FALSE),
                                "</br><b>", paste0(comp1.name, " ", df.vars$res1.sig.col), ":</b> ",
                                format(round(full.df$sig.x, 4), nsmall = 4),
                                "</br><b>", paste0(comp2.name, " ", df.vars$res2.sig.col), ":</b> ",
                                format(round(full.df$sig.y, 4), nsmall = 4),
                                "</br><b>", paste0(comp1.name, " ", df.vars$res1.expr.col),":</b> ",
                                format(round(full.df$exp.x, 2), nsmall = 2),
                                "</br><b>", paste0(comp2.name, " ", df.vars$res2.expr.col),":</b> ",
                                format(round(full.df$exp.y, 2), nsmall = 2))

  full.df <- as.data.frame(full.df)
  full.df <- full.df[order(full.df$order),]

  # Add plot border.
  ay <- list(
    showline = TRUE,
    mirror = TRUE,
    linecolor = toRGB("black"),
    linewidth = 0.5,
    title = paste0(comp2.name, "\n", df.vars$res2.lfc.col),
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
    title = paste0(comp1.name, "\n", df.vars$res1.lfc.col),
    range = list(-xlim, xlim),
    showgrid = FALSE,
    layer = "below traces",
    ticks = "outside",
    zerolinewidth = 0.5
  )

  fig <- plot_ly(full.df, x = ~lfc.x,
                 y = ~lfc.y,
                 customdata = ~Gene,
                 type = "scatter",
                 mode = "markers",
                 marker = list(color = ~col,
                               size = ~cex,
                               symbol = ~sh,
                               line = list(color = ~col),
                               opacity = opacity),
                 text = ~hover.string,
                 hoverinfo = "text",
                 source = source) %>%
    config(edits = list(annotationPosition = TRUE,
                        annotationTail = TRUE),
           toImageButtonOptions = list(format = "svg"),
           displaylogo = FALSE)

  if (regr) {
    fig <- fig %>%
      add_annotations(
        x= 0,
        y= 1,
        xref = "paper",
        yref = "paper",
        text = regr.anno,
        showarrow = FALSE,
        font = list(size = 10)
      )
  }

  if (!is.null(genes.labeled)) {
    fig <- fig %>%
      layout(xaxis = ax,
             yaxis = ay,
             showlegend = FALSE, shapes = list(regr.line)) %>%
      add_annotations(x = genes.labeled$x, y = genes.labeled$y, text = genes.labeled$customdata,
                      font = list(size = label.size, family = "Arial"), arrowside = "none")
  } else {
    fig <- fig %>% layout(xaxis = ax,
                   yaxis = ay, showlegend = FALSE, shapes = list(regr.line))
  }

  if (webgl) {
    fig <- fig %>% toWebGL()
  }

  fig
}


# Used to get column names that may differ between results.
.get_plot_vars <- function(res1, res2, sig.col, lfc.col, expr.col) {
  res1.names <- colnames(res1[[1]])
  res2.names <- colnames(res2[[1]])

  out <- list()
  out$res1.sig.col <- sig.col
  out$res2.sig.col <- sig.col
  out$res1.lfc.col <- lfc.col
  out$res2.lfc.col <- lfc.col
  out$res1.expr.col <- expr.col
  out$res2.expr.col <- expr.col

  # If column names are provided, assume they're the same for all results dataframes.
  if (is.null(sig.col)) {
    if (!any(res1.names %in% c("padj", "FDR", "svalue", "adj.P.Val"))) {
      stop("Cannot determine significance column, please provide the column name to sig.col")
    } else {
      out$res1.sig.col <- res1.names[res1.names %in% c("padj", "FDR", "svalue", "adj.P.Val")]
    }

    if (!any(res2.names %in% c("padj", "FDR", "svalue", "adj.P.Val"))) {
      stop("Cannot determine significance column, please provide the column name to sig.col")
    } else {
      out$res2.sig.col <- res2.names[res2.names %in% c("padj", "FDR", "svalue", "adj.P.Val")]
    }
  }

  if (is.null(lfc.col)) {
    if (!any(res1.names %in% c("log2FoldChange", "logFC", "LFC"))) {
      stop("Cannot determine fold change column, please provide the column name to lfc.col")
    } else {
      out$res1.lfc.col <- res1.names[res1.names %in% c("log2FoldChange", "logFC", "LFC")]
    }

    if (!any(res2.names %in% c("log2FoldChange", "logFC", "LFC"))) {
      stop("Cannot determine fold change column, please provide the column name to lfc.col")
    } else {
      out$res2.lfc.col <- res2.names[res2.names %in% c("log2FoldChange", "logFC", "LFC")]
    }
  }

  if (is.null(expr.col)) {
    if (!any(res1.names %in% c("baseMean", "logCPM", "AveExpr"))) {
      message("Cannot determine average expression column")
    } else {
      out$res1.expr.col <- res1.names[res1.names %in% c("baseMean", "logCPM", "AveExpr")]
    }

    if (!any(res2.names %in% c("baseMean", "logCPM", "AveExpr"))) {
      message("Cannot determine average expression column")
    } else {
      out$res2.expr.col <- res2.names[res2.names %in% c("baseMean", "logCPM", "AveExpr")]
    }
  }

  out
}
