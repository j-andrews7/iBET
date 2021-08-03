# Make the heatmap for differentially expressed genes under certain cutoffs.
.make_heatmap <- function(mat, res, anno, bm_col_func, lgc_col_func, fdr = 0.05, base_mean = 0, log2fc = 0, row_km = 0) {

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
            name = "log10(baseMean+1)", col = baseMean_col_fun, show_column_names = FALSE) +
    Heatmap(res$log2FoldChange[l], show_row_names = FALSE, width = unit(5, "mm"),
            name = "log2FoldChange", col = log2fc_col_fun, show_column_names = FALSE)
  ht <- draw(ht, merge_legend = TRUE)
  ht
}

.make_maplot <- function(res, ylim, highlight = NULL) {

  col <- rep("#00000020", nrow(res))
  cex <- rep(0.5, nrow(res))
  names(col) <- rownames(res)
  names(cex) <- rownames(res)
  if(!is.null(highlight)) {
    col[highlight] = "red"
    cex[highlight] = 1
  }
  x <- res$baseMean
  y <- res$log2FoldChange
  y[y > ylim] <- ylim
  y[y < -ylim] <- -ylim
  col[col == "red" & y < 0] <- "darkgreen"
  par(mar = c(4, 4, 1, 1))

  suppressWarnings(
    plot(x, y, col = col,
         pch = ifelse(res$log2FoldChange > ylim, 2, ifelse(res$log2FoldChange < -ylim, 6, 16)),
         cex = cex, log = "x",
         xlab = "baseMean", ylab = "log2 fold change")
  )
}

# make the volcano plot with some genes highlighted
.make_volcano <- function(res, highlight = NULL) {

  # Adjust for potential differences in the results table.
  sig.term <- "padj"
  if("svalue" %in% colnames(res)) {
    sig.term <- "svalue"
  }

  max.lfc <- max(abs(res$log2FoldChange)) + 0.2
  xlim <- c(-max.lfc, max.lfc)
  col <- rep("#00000020", nrow(res))
  cex <- rep(0.5, nrow(res))
  names(col) <- rownames(res)
  names(cex) <- rownames(res)
  if(!is.null(highlight)) {
    col[highlight] <- "red"
    cex[highlight] <- 1
  }
  x <- res$log2FoldChange
  y <- -log10(res[[sig.term]])
  ylim <- c(0, max(y[!is.infinite(y) & !is.na(y)]))
  y[is.infinite(y)] <- max(y[!is.infinite(y) & !is.na(y)]) + 1
  col[col == "red" & x < 0] <- "darkgreen"
  par(mar = c(4, 4, 1, 1))

  suppressWarnings(
    plot(x, y, col = col,
         pch = ifelse(y > max(ylim), 2, 16),
         cex = cex,
         xlab = "log2 fold change", ylab = paste0("-log10(", sig.term,")"),
         xlim = xlim,
         ylim = ylim)
  )
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
