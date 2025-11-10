# iBET - interactive Bioinformatics Exploratory Tools

Generating static reports from R markdown documents is a common way to distribute bioinformatics analyses with additional context and commentary. However, interactivity within these reports is limited. While one can embed minimally interactive figures in HTML reports, there is no way to alter the underlying data or plot by the end user.

[Shiny](https://shiny.rstudio.com/) provides an avenue to overcome this limitation, as [Shiny widgets](https://bookdown.org/yihui/rmarkdown/shiny-widgets.html) can be designed to interactively perform or display the results of an analysis. This comes with the downside that the Rmd file and necessary code must be hosted on a server that can run Shiny applications and interactive documents. [Options](https://www.rstudio.com/products/shiny/shiny-server/) include RStudio Connect, [shiny.apps.io](shiny.apps.io) (which has a free tier), and Shiny Server.

[iBET](https://github.com/j-andrews7/iBET) (interactive Bioinformatics Exploratory Tools) is an R package that contains drop-in Shiny widgets that run and/or display the results of common analyses performed with bioinformatics R packages. These empower end-users to fully immerse themselves in the data through interactive alterations of plot and analysis parameters. They are excellent tools for collaborative data wading and exploration and function well both as stand-alone apps or within Rmd files hosted on Shiny Servers.

## Installation

This package is in development - it may break at any time and contain unstable or untested features. A stable version will be submitted to CRAN or Bioconductor once the initially planned features are completed.

To install the package via Github:

```
install.packages("devtools")
devtools::install_github("j-andrews7/iBET")
```

## Usage

**iBET** currently contains three interactive widgets - `shinyPCAtools`, `shinyDE`, and `shinyDECorr`. They can be dropped into Rmd documents or ran directly within RStudio as shown below.

I **highly** recommend altering the width of your Rmd report by using a CSS block at the top of your document (right after the YAML header). This will use much more of the page, which makes using the widgets much easier. The following works well on most screens with no scaling:

```
```{css, echo=FALSE}
body {
  max-width: 1850px !important;
}
div.main-container {
  max-width: 1850px !important;
  width: 1850px !important;
  margin-left: auto !important;
  margin-right: auto !important;
}
.toc-content {
  padding-right: 0px;
  padding-left: 0px;
  margin-left: 300px;
  max-width: 1550px !important;
  width: 1550px !important;
}

```

### Load Example Data

```
shh <- suppressPackageStartupMessages
shh(library("airway"))
shh(library("magrittr"))
shh(library("DESeq2"))
shh(library("PCAtools"))
shh(library("iBET"))

data("airway")
airway$dex <- relevel(airway$dex, ref = "untrt")
dds <- DESeqDataSet(airway, design = ~ cell)
dds <- DESeq(dds)
vst <- assay(vst(dds))
```

### Interactive PCA via `PCAtools`

[PCAtools](https://bioconductor.org/packages/release/bioc/html/PCAtools.html) is a straight-forward package for principle component analysis.

```
shinyPCAtools(vst, metadata = colData(dds), annot.by = c("cell", "dex"), 
              color.by = "dex", shape.by = "cell", scale = TRUE, height = 850)
```


#### Larger (scRNA) Dataset

`shinyPCAtools` Works well on larger datasets, though it takes significant time for the PCA to run. Removing a much greater percentage of features is highly recommended. 

```
shh(library("scRNAseq"))
shh(library("scran"))
shh(library("scuttle"))
shh(library("scater"))
data <- ZeiselBrainData()
 
# # Remove genes expressed in few cells.
data <- data[rowMeans(counts(data) != 0) > 0.05, ]
data <- computeSumFactors(data, cluster = quickCluster(data))
data <- logNormCounts(data)
data$Cell.Type <- factor(data$level1class)

shinyPCAtools(logcounts(data), metadata = colData(data), annot.by = c("Cell.Type", "tissue", "age"), 
              color.by = "Cell.Type", removeVar = 0.9, scale = TRUE)
```

## Interactive Differential Expression Visualization

The `shinyDE` function provides an interactive widget for exploring differential expression results from any analysis tool (DESeq2, edgeR, limma, etc.). It was heavily inspired by the `interactivate` function from the [InteractiveComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/InteractiveComplexHeatmap.html) package.

Gene labels can be added to the MAplot and volcano plot by clicking a point. The labels can also be dragged around, though adding labels will reset the position, so it's recommended to add all labels prior to re-positioning them. Gene sets can be highlighted easily if provided.

Multiple comparisons can also be provided for switching between analyses quickly.

### Using with DESeq2 Results

```{r de, message = FALSE, warning = FALSE}
library(DESeq2)
deseq.res1 <- results(dds, contrast = c("dex", "trt", "untrt"))
dds <- DESeqDataSet(airway, design = ~ cell)
dds <- DESeq(dds)
deseq.res2 <- results(dds, contrast = c("cell", "N080611", "N052611"))
deseq.res3 <- results(dds, contrast = c("cell", "N61311", "N080611"))
deseq.res4 <- results(dds, contrast = c("cell", "N080611", "N61311"))
res <- list("trt v untrt" = as.data.frame(deseq.res1), 
            "N080611vN052611" = as.data.frame(deseq.res2), 
            "N61311vN080611" = as.data.frame(deseq.res3), 
            "N080611vN61311" = as.data.frame(deseq.res4))

# Get expression matrix and metadata
mat <- assay(vst(dds))
metadata <- as.data.frame(colData(dds))

shinyDE(mat = mat, res = res, metadata = metadata, 
        genesets = hs.msig, annot.by = c("cell", "dex"))
```

### Using with edgeR Results

```{r de_edger, message = FALSE, warning = FALSE, eval = FALSE}
library(edgeR)
y <- DGEList(counts = counts(dds))
y <- calcNormFactors(y)
design <- model.matrix(~ dex, data = colData(dds))
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
res_edger <- topTags(qlf, n = Inf)$table

# Use with shinyDE
shinyDE(
  mat = cpm(y, log = TRUE), 
  res = res_edger,
  metadata = colData(dds),
  lfc.col = "logFC",
  abundance.col = "logCPM", 
  sig.col = "FDR",
  annot.by = c("dex")
)
```

## Correlating DE Results

This widget creates scatter plots for arbitrary combinations of DE results. It takes a named list of dataframes as input and will plot scatter plots of log2 Fold Change values for each gene. It accepts up to 4 results dataframes and will generate plots for all combinations of them along with a regression line and correlation testing. Points are colored by significance in each DE analysis. Results dataframes from edgeR, limma, and DESeq2 will automatically have the appropriate columns used for plotting, but users can also provide the column names for their significance value or fold change column as necessary.

Gene labels can be added to a plot by clicking a point. The labels can also be dragged around, though adding labels to a plot will reset the label positions for said plot, so it's recommended to add all labels prior to re-positioning them.

```{r}
shinyDECorr(res, genesets = hs.msig, height = 800)
```

