# iBET - Interactive Bioconductor Exploratory Tools

Generating static reports from R markdown documents is a common way to distribute bioinformatics analyses with additional context and commentary. However, interactivity within these reports is limited. While one can embed minimally interactive figures in HTML reports, there is no way to alter the underlying data or plot by the end user.

[Shiny](https://shiny.rstudio.com/) provides an avenue to overcome this limitation, as [Shiny widgets](https://bookdown.org/yihui/rmarkdown/shiny-widgets.html) can be designed to interactively perform or display the results of an analysis. This comes with the downside that the Rmd file and necessary code must be hosted on a server that can run Shiny applications and interactive documents. [Options](https://www.rstudio.com/products/shiny/shiny-server/) include RStudio Connect, [shiny.apps.io](shiny.apps.io) (which has a free tier), and Shiny Server.

[iBET](https://github.com/j-andrews7/iBET) (interactive Bioconductor Exploratory Tools) is an R package that contains drop-in Shiny widgets that run and/or display the results of common analyses performed with Bioconductor R packages. These empower end-users to fully immerse themselves in the data through interactive alterations of plot and analysis parameters. Hopefully, it helps to alleviate some of the frustration that comes with frequent back-and-forth interactions to generate and alter figures.

## Installation

This package is in development - it may break at any time and contain unstable or untested features. A stable version will be submitted to CRAN or Bioconductor once the initially planned features are completed.

To install the package via Github:

```
install.packages("devtools")
devtools::install_github("j-andrews7/iBET")
```

## Usage

**iBET** currently contains two interactive widgets - `shinyPCAtools` and `shinyDESeq2`. They can be dropped into Rmd documents or ran directly within RStudio as shown below.

### Load Example Data

```
shh <- suppressPackageStartupMessages
shh(library("airway"))
shh(library("magrittr"))
shh(library("org.Hs.eg.db"))
shh(library("DESeq2"))
shh(library("PCAtools"))
shh(library("iBET"))
shh(library("scRNAseq"))
shh(library("scran"))
shh(library("scuttle"))
shh(library("scater"))

data('airway')
airway$dex %<>% relevel('untrt')

# Swap to gene symbols.
ens <- rownames(airway)
symbols <- mapIds(org.Hs.eg.db, keys = ens, column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(airway), names(symbols))]
rownames(airway) <- symbols
keep <- !is.na(rownames(airway))
airway <- airway[keep,]

# vst counts and make DESeqDataSet.
dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds)
vst <- assay(vst(dds))
```

### Interactive PCA via `PCAtools`

[PCAtools](https://bioconductor.org/packages/release/bioc/html/PCAtools.html) is a straight-forward package for principle component analysis.

```
shinyPCAtools(vst, metadata = colData(airway), annot.by = c("SampleName", "cell", "dex"), color.by = "dex", shape.by = "cell", height = 850)
```


#### Larger (scRNA) Dataset

`shinyPCAtools` Works well on larger datasets, though it takes significant time for the PCA to run. Removing a much greater percentage of features is highly recommended. 

```
data <- ZeiselBrainData()

# Remove genes expressed in few cells.
data <- data[rowMeans(counts(data) != 0) > 0.05, ]
data <- computeSumFactors(data, cluster = quickCluster(data))
data <- logNormCounts(data)
data$Cell.Type <- factor(data$level1class)

shinyPCAtools(logcounts(data), metadata = colData(data), annot.by = c("Cell.Type", "tissue", "age"), color.by = "Cell.Type", removeVar = 0.9)
```

### Interactive Differential Expression via `DESeq2`

This widget was heavily inspired by the `interactivate` function from the [InteractiveComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/InteractiveComplexHeatmap.html) package. It wraps [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and will run differential expression analysis if not provided a results dataframe as well.

```
shinyDESeq2(dds, coef = "dex_trt_vs_untrt", annot.by = c("cell", "dex"))
```


