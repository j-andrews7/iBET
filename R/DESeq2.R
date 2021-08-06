#' Create an interactive Shiny app for DESeq2 results exploration
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @import DESeq2
#' @import InteractiveComplexHeatmap
#' @import ComplexHeatmap
#' @import DT
#' @importFrom GetoptLong qq
#' @import circlize
#' @importFrom grid unit grid.newpage grid.text
#' @importFrom SummarizedExperiment colData assay
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics par
#' @importFrom stats quantile
#' @importFrom plotly ggplotly plotlyOutput renderPlotly toWebGL
#' @import ggplot2
#'
#' @param dds A \code{\link[DESeq2]{DESeqDataSet}} object.
#' @param res The object returned by \code{\link[DESeq2]{results}} or \code{\link[DESeq2]{lfcShrink}} (recommended).
#'   If not provided, it will be generated from the \code{dds} object via \code{\link[DESeq2]{lfcShrink}} and \code{coef}.
#' @param coef A string indicating the coefficient name for which results will be generated.
#'   If not provided and \code{res} is \code{NULL}, the first non-intercept coefficient
#'   provided by \code{\link[DESeq2]{resultsNames}} will be used.
#' @param annot.by A string or character vector containing the names of sample metadata variables to be used as heatmap annotations.
#'   If not provided, all variables used in the \code{design} will be used.
#' @param use.lfcShrink Boolean indicating whether \code{\link[DESeq2]{lfcShrink}} should be used during results
#'   generation. These are useful for ranking and visualization without the need for arbitrary filtering of low count genes.
#' @param lfcThreshold Number passed to \code{\link[DESeq2]{lfcShrink}} or \code{\link[DESeq2]{results}} to set an LFC threshold.
#'   p-values, s-values, and adjusted p-values returned will be for whether the LFC is significantly greater in absolute value than the threshold.
#'   Ignored if \code{res} is not \code{NULL}.
#' @param use.vst Boolean indicating whether \code{\link[DESeq2]{vst}} transformed counts should be used for heatmap generation (recommended). If \code{FALSE},
#'   normalized counts will be used. The variance stabilizing transformation helps to reduce increased variance seen for low counts in log space.
#'   This generally results in a better looking heatmap with fewer outliers.
#' @param h.id String indicating unique ID for interactive heatmaps.
#'   Required if multiple apps are run within the same Rmd file.
#'
#' @return A Shiny app containing an InteractiveComplexHeatmap, MAplot, and volcano plot that are interconnected.
#'
#'
#' @seealso
#' \code{\link[DESeq2]{results}}, \code{\link[DESeq2]{lfcShrink}}, \code{\link[DESeq2]{resultsNames}}.
#'
#' @author Jared Andrews, based heavily on code by Zuguang Gu.
#' @export
shinyDESeq2 <- function(dds, res = NULL, coef = NULL, annot.by = NULL,
                        use.lfcShrink = TRUE, lfcThreshold = 0, use.vst = TRUE, h.id = "ht1") {

  # If gene dispersions not yet calculated, calculate them.
  if(is.null(body(dds@dispersionFunction))) {
    dds <- DESeq2::DESeq(dds)
  }

  # Get contrast name if results are not provided.
  if (is.null(res)) {
    if(is.null(coef)) {
      coef <- resultsNames(dds)[2]
    }

    if (use.lfcShrink) {
      res <- DESeq2::lfcShrink(dds, coef = coef, lfcThreshold = lfcThreshold)
    } else {
      res <- DESeq2::results(dds, name = coef, lfcThreshold = lfcThreshold)
    }
  }

  if (use.vst) {
    mat <- SummarizedExperiment::assay(DESeq2::vst(dds))
  } else{
    mat <- as.matrix(DESeq2::counts(dds, normalized = TRUE))
  }

  anno <- SummarizedExperiment::colData(dds)

  # Change to use `annot.by`.
  if (!is.null(annot.by)) {
    anno <- anno[, annot.by, drop = FALSE]
  } else {
    anno <- anno[, all.vars(dds@design), drop = FALSE]
  }

  l <- sapply(anno, function(x) (is.factor(x) || is.character(x)) && !any(duplicated(x))) | sapply(anno, function(x) length(unique(x)) == 1)
  anno <- anno[, !l, drop = FALSE]

  env <- new.env()

  qa <- quantile(log10(res$baseMean + 1), 0.99)
  baseMean_col_fun <- circlize::colorRamp2(c(0, qa/2, qa), c("blue", "white", "red"))

  qa <- quantile(abs(abs(res$log2FoldChange[!is.na(res$log2FoldChange)])), 0.99)
  log2fc_col_fun <- circlize::colorRamp2(c(-qa, 0, qa), c("green", "white", "red"))

  environment(.make_heatmap) <- env

  # A self-defined action to respond brush event. It updates the MA-plot, the volcano plot
  # and a table which contains DESeq2 results for the selected genes.
  .brush_action <- function(df, input, output, session) {

    row_index <- unique(unlist(df$row_index))
    selected <- env$row_index[row_index]

    output[["ma_plot"]] <- renderPlotly({
      .make_maplot(res = res, ylim = input$ma.y, highlight = selected)
    })

    output[["volcano_plot"]] <- renderPlotly({
      .make_volcano(res = res,xlim = input$vol.x, ylim = input$vol.y, highlight = selected)
    })

    output[["res_table"]] <- DT::renderDataTable({
      # Adjust output table columns based on results table.
      if ("svalue" %in% colnames(res)) {
        third <- "svalue"
      } else {
        third <- "padj"
      }

      DT::formatRound(DT::datatable(as.data.frame(res[selected, c("baseMean", "log2FoldChange", third)]), rownames = TRUE), columns = 1:3, digits = 3)
    })
  }

  .click_action <- function(df, input, output, session) {
    row_index <- unique(unlist(df$row_index))
    selected <- env$row_index[row_index]

    output[["ma_plot"]] <- renderPlotly({
      .make_maplot(res = res, ylim = input$ma.y, highlight = selected)
    })

    output[["volcano_plot"]] <- renderPlotly({
      .make_volcano(res = res,xlim = input$vol.x, ylim = input$vol.y, highlight = selected)
    })

    output[["res_table"]] <- DT::renderDataTable({
      # Adjust output table columns based on results table.
      if ("svalue" %in% colnames(res)) {
        third <- "svalue"
      } else {
        third <- "padj"
      }

      DT::formatRound(DT::datatable(as.data.frame(res[selected, c("baseMean", "log2FoldChange", third)]),
                                    rownames = TRUE), columns = 1:3, digits = 5)
    })
  }

  # The dashboard body contains three columns:
  # 1. the original heatmap
  # 2. the sub-heatmap and the default output
  # 3. the self-defined output
  body <- mainPanel(width = 10,
    fluidRow(
      column(width = 4,
        h3("Differential heatmap"),
        originalHeatmapOutput(h.id, height = 600, width = 300, containment = TRUE)
      ),
      column(width = 4,
        id = "column2",
        h3("Sub-heatmap"),
        subHeatmapOutput(h.id, title = NULL, width = 300, containment = TRUE),
        h3(title = "Output"),
        HeatmapInfoOutput(h.id, title = NULL, width = 300),
        h3("Result table of the selected genes"),
        DT::dataTableOutput("res_table")
      ),
      column(width = 4,
        h3("MA-plot"),
        plotlyOutput("ma_plot", height = "350px", width = "300px"),
        htmlOutput("ma_plot_selected"),
        h3("Volcano plot"),
        plotlyOutput("volcano_plot", height = "350px", width = "300px"),
        htmlOutput("volcano_plot_selected"),
      )
    )
  )

  # Side bar contains settings for certain cutoffs to select significant genes.
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(width = 2,
        tags$label(HTML(qq("Comparison: <code style='font-weight:normal;'>@{paste(coef, collapse = ' ')}</code>")), class = "shiny-input-container", style = "font-size:1.2em;"),
        hr(style="margin:2px;"),
        numericInput("fdr", label = "Significance threshold:", value = 0.05, step = 0.001, min = 0.0001),
        numericInput("base_mean", label = "Minimal base mean:", value = 0, step = 1),
        numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 0, step = 0.1, min = 0),
        numericInput("km", label = "Number of k-means groups. Set to 0 to suppress k-means clustering:", value = 2, step = 1),
        numericInput("ma.y", label = "MAplot y-axis limits:", value = 3, step = 0.1, min = 0.1),
        numericInput("vol.x", label = "Volcano plot x-axis limits:", value = 3, step = 0.1, min = 0.1),
        numericInput("vol.y", label = "Volcano plot y-axis limits:", value = 50, step = 0.5, min = 1),
        actionButton("filter", label = "Generate heatmap")
      ),
      body
    )
  )

  # makeInteractiveComplexHeatmap() is put inside observeEvent() so that changes on the cutoffs can regenerate the heatmap.
  server <- function(input, output, session) {
    observeEvent(input$filter, {
      pdf(NULL)
      ht <- .make_heatmap(mat, res, anno, baseMean_col_fun, log2fc_col_fun,
                          fdr = as.numeric(input$fdr), base_mean = input$base_mean, log2fc = input$log2fc,
                        row_km = input$km)
      dev.off()

      if(!is.null(ht)) {
        makeInteractiveComplexHeatmap(input, output, session, ht, h.id,
                                      brush_action = .brush_action, click_action = .click_action)
      } else {
        # The ID for the heatmap plot is encoded as @{heatmap_id}_heatmap, thus, it is ht_heatmap here.
        output[[paste0(h.id, "_heatmap")]] <- renderPlot({
          grid.newpage()
          grid.text("No row exists after filtering.")
        })
      }
    }, ignoreNULL = FALSE)

    # observeEvent(c(input$ma.y), {
    #   #selected <- env$row_index[row_index]
    #   map <- debounce(.make_maplot(res = res, ylim = input$ma.y, highlight = NULL), 2000)
    #   output[["ma_plot"]] <- renderPlotly({
    #     map
    #   })
    # }, ignoreInit = TRUE)
  }

  shinyApp(ui, server)
}
