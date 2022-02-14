# Generate easier columns for plotting for various data summaries.
.gene_ingress <- function(df, sig.thresh, lfc.thresh, essential.genes = NULL, depmap.genes = NULL) {
  
  if (!is.null(essential.genes)) {
    df$essential <- df$id %in% essential.genes
  } 
  
  if (!is.null(depmap.genes)) {
    df$DepMap_CRISPR_Essential <- df$id %in% depmap.genes$Gene[depmap.genes$Dataset == "DependencyEnum.Chronos_Combined" & 
                                                                 depmap.genes$Common.Essential == "True"]
    df$DepMap_CRISPR_Selective <- df$id %in% depmap.genes$Gene[depmap.genes$Dataset == "DependencyEnum.Chronos_Combined" & 
                                                                 depmap.genes$Strongly.Selective == "True"]
    
    df$DepMap_RNAi_Essential <- df$id %in% depmap.genes$Gene[depmap.genes$Dataset == "DependencyEnum.RNAi_merged" & 
                                                               depmap.genes$Common.Essential == "True"]
    df$DepMap_RNAi_Selective <- df$id %in% depmap.genes$Gene[depmap.genes$Dataset == "DependencyEnum.RNAi_merged" & 
                                                               depmap.genes$Strongly.Selective == "True"]
  }
  
  df$LFC <- as.numeric(df$`neg|lfc`)
  
  df$RRAscore <- apply(df, 1, function(x) {
    x <- split(unname(x),names(x))
    as.numeric(ifelse(as.numeric(x$`neg|score`) < as.numeric(x$`pos|score`), x$`neg|score`, x$`pos|score`))
  })
  
  df$FDR <- apply(df, 1, function(x) {
    x <- split(unname(x),names(x))
    as.numeric(ifelse(as.numeric(x$`neg|fdr`) < as.numeric(x$`pos|fdr`), x$`neg|fdr`, x$`pos|fdr`))
  })
  
  df$hit_type <- apply(df, 1, function(x) {
    x <- split(unname(x),names(x))
    ifelse(as.numeric(x$`neg|fdr`) < sig.thresh & as.numeric(x$LFC) < -as.numeric(lfc.thresh), "neg", 
           ifelse(as.numeric(x$`pos|fdr`) < sig.thresh & as.numeric(x$LFC) > as.numeric(lfc.thresh), "pos", NA))
  })
  
  df$pval <- apply(df, 1, function(x) {
    x <- split(unname(x),names(x))
    as.numeric(ifelse(as.numeric(x$`neg|p-value`) < as.numeric(x$`pos|p-value`), x$`neg|p-value`, x$`pos|p-value`))
  })
  
  df$goodsgrna <- apply(df, 1, function(x) {
    x <- split(unname(x),names(x))
    as.integer(ifelse(as.numeric(x$`neg|score`) < as.numeric(x$`pos|score`), x$`neg|goodsgrna`, x$`pos|goodsgrna`))
  })
  
  df$Rank <- rank(df$LFC)
  
  df$RandomIndex <- sample(1:nrow(df), nrow(df))

  df
}