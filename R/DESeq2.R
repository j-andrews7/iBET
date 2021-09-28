#' Create an interactive Shiny app for DESeq2 analysis and results exploration
#'
#' This shiny app is composed of three tabs - an interactive heatmap, MA and volcano plots, and a table of full
#' differential expression results. The interactive heatmap will generate a sub-heatmap for selected rows/columns.
#'
#' Gene labels can be added to the MAplot and volcano plot by clicking a point. The labels can also be dragged around,
#' though adding labels will reset the positions, so it's recommended to add all labels prior to re-positioning them.
#'
#' @details Note that significance values of 0 will always be pushed above the y-limit of the volcano plot, as they are infinite values after log transformation.
#'
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @import DESeq2
#' @import InteractiveComplexHeatmap
#' @import ComplexHeatmap
#' @import DT
#' @importFrom GetoptLong qq
#' @import circlize
#' @importFrom grid unit grid.newpage grid.text gpar
#' @importFrom SummarizedExperiment colData assay
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics par
#' @importFrom stats quantile
#' @importFrom plotly ggplotly plotlyOutput renderPlotly toWebGL plot_ly layout add_annotations config toRGB event_data
#' @import ggplot2
#' @import shinydashboard
#' @import dashboardthemes
#' @importFrom shinyWidgets prettyCheckbox dropdownButton tooltipOptions
#' @importFrom shinycustomloader withLoader
#' @importFrom shinyjqui jqui_resizable
#' @importFrom colourpicker colourInput
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
#' @param samples.use Vector of Boolean values indicating which samples should be included in the visualizations.
#'   Must be the same length as the number of samples in the \code{dds} object.
#' @param h.id String indicating unique ID for interactive heatmaps.
#'   Required if multiple apps are run within the same Rmd file.
#' @param height Number indicating height of app in pixels.
#'
#' @return A Shiny app containing interconnected InteractiveComplexHeatmap, MAplot, and volcano plots along with full DE results.
#'
#'
#' @seealso
#' \code{\link[DESeq2]{results}}, \code{\link[DESeq2]{lfcShrink}}, \code{\link[DESeq2]{resultsNames}}.
#'
#' @author Jared Andrews, based on code by Zuguang Gu.
#' @export
shinyDESeq2 <- function(dds, res = NULL, coef = NULL, annot.by = NULL,
                        use.lfcShrink = TRUE, lfcThreshold = 0, use.vst = TRUE, samples.use = NULL,
                        h.id = "ht1", height = 800) {

  # If gene dispersions not yet calculated, calculate them.
  if(is.null(body(dds@dispersionFunction))) {
    dds <- DESeq2::DESeq(dds)
  }

  # Get contrast name if results are not provided.
  if (is.null(res)) {
    if(is.null(coef)) {
      coef <- resultsNames(dds)[resultsNames(dds) != "Intercept"][1]
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

  # Filter samples.
  if (!is.null(samples.use)) {
    mat <- mat[,samples.use]
    anno <- anno[samples.use,]
  }

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
  log2fc_col_fun <- circlize::colorRamp2(c(-qa, 0, qa), c("blue", "white", "red"))

  if ("svalue" %in% colnames(res)) {
    sig.term <- "svalue"
  } else {
    sig.term <- "padj"
  }

  environment(.make_heatmap) <- env

  # A self-defined action to respond brush event. It updates the MA-plot, the volcano plot
  # and a table which contains DESeq2 results for the selected genes.
  .brush_action <- function(df, input, output, session) {

    row_index <- unique(unlist(df$row_index))
    selected <- env$row_index[row_index]

    output[["res_table"]] <- DT::renderDataTable({
      # Adjust output table columns based on results table.
      if ("svalue" %in% colnames(res)) {
        third <- "svalue"
      } else {
        third <- "padj"
      }

      DT::formatRound(DT::datatable(as.data.frame(res[selected, c("baseMean", "log2FoldChange", sig.term)]),
                                    rownames = TRUE, options = list(lengthMenu = c(5, 10, 25),
                                                                    pageLength = 20)),
                      columns = 1:3, digits = 5) %>%
        DT::formatStyle(0, target = "row", lineHeight = '40%')
    })
  }

  .click_action <- function(df, input, output, session) {
    row_index <- unique(unlist(df$row_index))
    selected <- env$row_index[row_index]

    output[["res_table"]] <- DT::renderDataTable({
      # Adjust output table columns based on results table.
      if ("svalue" %in% colnames(res)) {
        third <- "svalue"
      } else {
        third <- "padj"
      }

      DT::formatRound(DT::datatable(as.data.frame(res[selected, c("baseMean", "log2FoldChange", sig.term)]),
                                    rownames = TRUE, options = list(lengthMenu = c(5, 10, 25),
                                                                  pageLength = 20)),
                      columns = 1:3, digits = 5) %>%
        DT::formatStyle(0, target = "row", lineHeight = '40%')
    })
  }


  body <- mainPanel(width = 10,
    tabsetPanel(
      tabPanel("Heatmap",
        fluidRow(
          column(width = 4,
            h3("Differential heatmap"),
            originalHeatmapOutput(h.id, height = 500, width = 400, containment = TRUE)
          ),
          column(width = 4,
            id = "column2",
            h3("Sub-heatmap"),
            subHeatmapOutput(h.id, title = NULL, width = 400, containment = TRUE),
            h3(title = "Output"),
            HeatmapInfoOutput(h.id, title = NULL, width = 400),
          ),
          column(width = 4,
                 h3("Results for the Selected Genes"),
                 div(DT::dataTableOutput("res_table"), style = "font-size:80%")
          )
        )
      ),
      tabPanel("MA & Volcano Plots",
        fluidRow(
          column(width = 6,
            br(),
            dropdownButton(
              tags$h3("Plot Settings"),
              colourInput("ma.down.color", "Down-regulated genes colour", value = "#0026ff"),
              colourInput("ma.up.color", "Up-regulated genes colour", value = "red"),
              colourInput("ma.insig.color", "Unchanged genes colour", value = "black"),
              numericInput("ma.y", label = "MAplot y-axis limits:", value = 5, step = 0.1, min = 0.1),
              prettyCheckbox("ma.fcline", label = "Show MAplot FC Threshold", value = TRUE,
                            animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
              circle = FALSE, label = strong("MA-Plot"), status = "danger", size = "lg", icon = icon("gear"),
              width = "300px", tooltip = tooltipOptions(title = "Click to change plot settings")
            ),
            withLoader(
              jqui_resizable(
                plotlyOutput("ma_plot", height = "550px", width = "550px")
              ),
              type = "html", loader = "dnaspin"
            )
          ),
          column(width = 6,
            br(),
            dropdownButton(
              tags$h3("Plot Settings"),
              colourInput("vol.down.color", "Down-regulated genes colour", value = "#0026ff"),
              colourInput("vol.up.color", "Up-regulated genes colour", value = "red"),
              colourInput("vol.insig.color", "Unchanged genes colour", value = "#A6A6A6"),
              numericInput("vol.x", label = "Volcano plot x-axis limits:", value = 5, step = 0.1, min = 0.1),
              numericInput("vol.y", label = "Volcano plot y-axis limits:", value = max(-log10(res[[sig.term]])),
                           step = 0.5, min = 1),
              prettyCheckbox("vol.fcline", label = "Show Volcano FC Threshold", value = TRUE,
                             animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
              prettyCheckbox("vol.sigline", label = "Show Volcano Signficance Threshold", value = TRUE,
                             animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
              circle = FALSE, label = strong("Volcano Plot"), status = "danger", size = "lg", icon = icon("gear"),
              width = "300px", tooltip = tooltipOptions(title = "Click to change plot settings")
            ),
            withLoader(jqui_resizable(plotlyOutput("volcano_plot", height = "550px", width = "550px")),
                       type = "html", loader = "dnaspin"),
          )
        )
      ),
      tabPanel("Full Results Table",
        div(DT::dataTableOutput("res_table_full"), style = "font-size:80%; margin:10px;")
      )
    )
  )

  # Side bar contains settings for certain cutoffs to select significant genes.
  ui <- dashboardPage(
    dashboardHeader(disable = TRUE),
    dashboardSidebar(disable = TRUE),
    dashboardBody(
      shinyDashboardThemes(
        theme = "onenote"
      ),
      sidebarLayout(
        sidebarPanel(width = 2,
          tags$label(HTML(qq("Comparison: <code style='font-weight:normal; font-size: 10px;'>@{paste(coef, collapse = ' ')}</code>")), class = "shiny-input-container", style = "font-size:1.2em;"),
          hr(style="margin:2px; background-color: #737373;"),
          numericInput("fdr", label = "Significance threshold:", value = 0.05, step = 0.001, min = 0.0001),
          numericInput("base_mean", label = "Minimal base mean:", value = 0, step = 1),
          numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 0, step = 0.1, min = 0),
          numericInput("row.km", label = "Number of row k-means groups:", value = 2, step = 1),
          numericInput("col.km", label = "Number of column k-means groups:", value = 0, step = 1),
          actionButton("filter", label = "Generate heatmap")
        ),
        body
      )
    )
  )

  # makeInteractiveComplexHeatmap() is put inside observeEvent() so that changes on the cutoffs can regenerate the heatmap.
  server <- function(input, output, session) {

    # Keep track of which genes have been clicked
    genes <- reactiveValues(ma = NULL, volc = NULL)

    # On click, the key field of the event data contains the gene symbol
    # Add that gene to the set of all "selected" genes
    observeEvent(event_data("plotly_click", source = paste0(h.id,"_ma")), {
      gene <- event_data("plotly_click", source = paste0(h.id,"_ma"))
      gene_old_new <- rbind(genes$ma, gene)
      keep <- gene_old_new[gene_old_new$customdata %in% names(which(table(gene_old_new$customdata)==1)),]

      if (nrow(keep) == 0) {
        genes$ma <- NULL
      } else {
        genes$ma <- keep
      }
    })

    observeEvent(event_data("plotly_click", source = paste0(h.id,"_volc")), {
      gene <- event_data("plotly_click", source = paste0(h.id,"_volc"))
      gene_old_new <- rbind(genes$volc, gene)
      keep <- gene_old_new[gene_old_new$customdata %in% names(which(table(gene_old_new$customdata)==1)),]

      if (nrow(keep) == 0) {
        genes$volc <- NULL
      } else {
        genes$volc <- keep
      }
    })

    # clear the set of genes when a double-click occurs
    observeEvent(event_data("plotly_doubleclick", source = paste0(h.id,"_ma")), {
      genes$ma <- NULL
    })

    observeEvent(event_data("plotly_doubleclick", source = paste0(h.id,"_volc")), {
      genes$volc <- NULL
    })

    output$ma_plot <- renderPlotly({
      req(genes)
      .make_maplot(res = res, ylim = input$ma.y, fc.thresh = input$log2fc,
                   fc.lines = input$ma.fcline, sig.thresh = input$fdr, h.id = h.id,
                   sig.term = sig.term, gs = genes$ma, up.color = input$ma.up.color,
                   down.color = input$ma.down.color, insig.color = input$ma.insig.color)
    })

    output$volcano_plot <- renderPlotly({
      req(genes)
      .make_volcano(res = res, xlim = input$vol.x, ylim = input$vol.y, fc.thresh = input$log2fc,
                    fc.lines = input$vol.fcline, sig.thresh = input$fdr, sig.line = input$vol.sigline,
                    h.id = h.id, sig.term = sig.term, gs = genes$volc, up.color = input$vol.up.color,
                    down.color = input$vol.down.color, insig.color = input$vol.insig.color)
    })

    output[["res_table_full"]] <- DT::renderDataTable({
      df <- as.data.frame(res)

      # Collect proper output columns..
      if (sig.term == "svalue") {
        cnames <- "svalue"
      } else if ("stat" %in% colnames(res)) {
        cnames <- c("padj", "pvalue", "stat")
      } else {
        cnames <- c("padj", "pvalue")
      }

      df <- df[, c("baseMean", "log2FoldChange", "lfcSE", cnames)]

      DT::formatRound(DT::datatable(df,
                                    rownames = TRUE,
                                    filter = "top",
                                    extensions = c("Buttons"),
                                    options = list(lengthMenu = c(5, 10, 25, 50),
                                                   pageLength = 25,
                                                   dom = 'Blfrtip',
                                                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))),
                      columns = c("baseMean", "log2FoldChange", cnames), digits = 5) %>%
        DT::formatStyle(1, target = "row", lineHeight = '40%')
    })

    observeEvent(input$filter, {
      pdf(NULL)
      ht <- .make_heatmap(mat, res, anno, baseMean_col_fun, log2fc_col_fun,
                          fdr = as.numeric(input$fdr), base_mean = input$base_mean, log2fc = input$log2fc,
                        row.km = input$row.km, col.km = input$col.km)
      dev.off()

      if (!is.null(ht)) {
        makeInteractiveComplexHeatmap(input, output, session, ht, h.id,
                                      brush_action = .brush_action, click_action = .click_action)
      } else {
        # The ID for the heatmap plot is encoded as @{heatmap_id}_heatmap.
        output[[paste0(h.id, "_heatmap")]] <- renderPlot({
          grid.newpage()
          grid.text("No row exists after filtering.")
        })
      }
    }, ignoreNULL = FALSE)

  }

  shinyApp(ui, server, options = list(height = height))
}
