# Make the heatmap for differentially expressed genes under certain cutoffs.
.make_heatmap <- function(mat, res, anno, bm_col_func, lfc_col_func, fdr = 0.05, base_mean = 0, log2fc = 0, row_km = 0) {

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
                show_row_names = FALSE, show_column_names = FALSE, row_km = row_km,
                column_title = paste0(sum(l), " significant genes with ", sig.term," < ", fdr),
                show_row_dend = FALSE) +
    Heatmap(log10(res$baseMean[l]+1), show_row_names = FALSE, width = unit(5, "mm"),
            name = "log10(baseMean+1)", col = bm_col_func, show_column_names = FALSE) +
    Heatmap(res$log2FoldChange[l], show_row_names = FALSE, width = unit(5, "mm"),
            name = "log2FoldChange", col = lfc_col_func, show_column_names = FALSE)
  ht <- draw(ht, merge_legend = TRUE)
  ht
}

.make_maplot <- function(res, ylim, highlight = NULL) {

  # Adjust for potential differences in the results table.
  sig.term <- "padj"
  if("svalue" %in% colnames(res)) {
    sig.term <- "svalue"
  }

  col <- rep("#00000020", nrow(res))
  cex <- rep(0.5, nrow(res))
  names(col) <- rownames(res)
  names(cex) <- rownames(res)
  if(!is.null(highlight)) {
    col[highlight] = "red"
    cex[highlight] = 1
  }
  res$x <- res$baseMean
  res$y <- res$log2FoldChange
  sh <- ifelse(res$log2FoldChange > ylim, 2, ifelse(res$log2FoldChange < -ylim, 6, 16))
  res$y[res$y > ylim] <- ylim
  res$y[res$y < -ylim] <- -ylim
  col[col == "red" & res$log2FoldChange < 0] <- "darkgreen"
  res$Gene <- rownames(res)

  res$hover.string <- paste("</br><b>Gene:</b> ", res$Gene,
                            "</br><b>log2 Fold Change:</b> ", format(round(res$log2FoldChange, 4), nsmall = 4),
                            "</br><b>", sig.term, ":</b> ", format(round(res[[sig.term]], 4), nsmall = 4),
                            "</br><b>baseMean (avg. Expression):</b> ", format(round(res$baseMean, 2), nsmall = 2))

  ggplotly(
    ggplot(as.data.frame(res), aes_string(x, y, text = hover.string)) + geom_point(col = col,
      shape = sh, size = cex) + xlab("baseMean") + ylab("log2 fold change") +
      ylim(-ylim, ylim) + scale_x_log10() + theme_bw(), tooltip = c("text")
  ) %>% toWebGL()
}

# make the volcano plot with some genes highlighted
.make_volcano <- function(res, xlim, ylim, highlight = NULL) {

  # Adjust for potential differences in the results table.
  sig.term <- "padj"
  if("svalue" %in% colnames(res)) {
    sig.term <- "svalue"
  }

  col <- rep("#00000020", nrow(res))
  cex <- rep(0.5, nrow(res))
  names(col) <- rownames(res)
  names(cex) <- rownames(res)
  if(!is.null(highlight)) {
    col[highlight] <- "red"
    cex[highlight] <- 1
  }
  res$x <- res$log2FoldChange
  res$y <- -log10(res[[sig.term]])
  #sh <- ifelse(res$y > ylim, 2, ifelse(res$x < -xlim, 60, ifelse(res$x > xlim, 62, 16)))
  sh <- ifelse(res$y > ylim, "triangle-up-open",
               ifelse(res$x < -xlim, "triangle-left-open",
                      ifelse(res$x > xlim, "triangle-right-open", 16)))

  res$y[res$y > ylim] <- ylim
  res$x[res$x > xlim] <- xlim
  res$x[res$x < -xlim] <- -xlim
  ylim <- c(0, ylim)

  col[col == "red" & res$x < 0] <- "darkgreen"

  res$Gene <- rownames(res)

  res$hover.string <- paste("</br><b>Gene:</b> ", res$Gene,
                            "</br><b>log2 Fold Change:</b> ", format(round(res$log2FoldChange, 4), nsmall = 4),
                            "</br><b>", sig.term, ":</b> ", format(round(res[[sig.term]], 6), nsmall = 6),
                            "</br><b>baseMean (avg. Expression):</b> ", format(round(res$baseMean, 2), nsmall = 2))

  ggplotly(
    ggplot(as.data.frame(res), aes_string(x, y, text = hover.string)) +
      geom_point(col = col, shape = sh, size = cex) +
      xlab("log2 fold change") + ylab(paste0("-log10(", sig.term,")")) +
      xlim(-xlim, xlim) + ylim(ylim) +
      theme_bw(), tooltip = c("text")
  ) %>% toWebGL()
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
