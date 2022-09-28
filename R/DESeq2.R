#' Create an interactive Shiny app for DESeq2 analysis and results exploration
#'
#' This shiny app is composed of three tabs - an interactive heatmap, MA and volcano plots, and a table of full
#' differential expression results. The interactive heatmap will generate a sub-heatmap for selected rows/columns.
#'
#' Gene labels can be added to the MAplot and volcano plot by clicking a point. The labels can also be dragged around,
#' though adding labels will reset the positions, so it's recommended to add all labels prior to re-positioning them.
#'
#' @details Note that significance values of 0 will always be pushed above the y-limit of the volcano plot,
#'   as they are infinite values after log transformation.
#'
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
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
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyjqui jqui_resizable
#' @importFrom shinyjs show useShinyjs hidden
#' @importFrom shinyBS tipify popify
#' @importFrom colourpicker colourInput
#' @importFrom methods is
#' @importFrom htmlwidgets saveWidget
#' @importFrom DESeq2 DESeq results lfcShrink resultsNames counts
#'
#' @param dds A \code{\link[DESeq2]{DESeqDataSet}} object.
#' @param res Either the object returned by \code{\link[DESeq2]{results}} or \code{\link[DESeq2]{lfcShrink}} (recommended)
#'   or a named list of such results. If a named list is provided, users will be
#'   able to choose between the provided results and \code{coef} will be ignored.
#'
#'   If not provided, it will be generated from the \code{dds} object via \code{\link[DESeq2]{lfcShrink}}
#'   or \code{\link[DESeq2]{results}} and \code{coef}.
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
#' @param genesets Optional named list containing genesets that can be interactively highlighted on the plots.
#'   The elements of the list should each be a geneset with gene identifiers matching those used in the results.
#' @param swap.rownames String. The column name of rowData(dds) \emph{and} res used to identify features instead of rownames(object).
#'
#'   Note if this column contains duplicates (e.g. gene symbols), they will be make unique via the addition of ".1", ".2", etc.
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
                        use.lfcShrink = TRUE, lfcThreshold = 0, use.vst = TRUE,
                        h.id = "ht1", genesets = NULL, swap.rownames = NULL, height = 800) {

  # Swap rownames if necessary.
  dds <- .swap_rownames(dds, swap.rownames = swap.rownames)

  # Check is res is individual or named list.
  multi.res <- FALSE
  res.list <- NULL
  if (is(res, "list")) {
    if (!is.null(names(res))) {
      multi.res <- TRUE
      # Swap rownames if necessary.
      res.list <- lapply(res, .swap_rownames,  swap.rownames = swap.rownames)
      res <- res.list[[1]]
    } else {
      stop("Results list elements should be named.")
    }
  }

  # Parameter validation.
  if (!is.null(genesets)) {
    if (is.null(names(genesets))) {
      stop("Genesets list must be named")
    } else if (!is(genesets, "list")) {
      stop("Genesets must be provided as a named list")
    }
  }

  if (!is.null(annot.by)) {
    if (!all(annot.by %in% names(SummarizedExperiment::colData(dds)))) {
      stop("Annotation variables not found in object colData")
    }
  }

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

  env <- new.env()

  qa <- quantile(log10(res$baseMean + 1), 0.99)
  baseMean_col_fun <- circlize::colorRamp2(c(0, qa/2, qa), c("blue", "white", "red"))

  qa <- quantile(abs(res$log2FoldChange[!is.na(res$log2FoldChange)]), 0.99)
  log2fc_col_fun <- circlize::colorRamp2(c(-qa, 0, qa), c("blue", "white", "red"))

  if ("svalue" %in% colnames(res)) {
    sig.term <- "svalue"
  } else {
    sig.term <- "padj"
  }

  environment(.make_heatmap) <- env

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
              fluidRow(
                column(width = 6,
                  tipify(colourInput("ma.down.color", "Down-reg colour", value = "#0026ff"),
                         "Color of down-regulated genes.", "right", options = list(container = "body")),
                  tipify(colourInput("ma.up.color", "Up-reg colour", value = "red"),
                         "Color of up-regulated genes.", "right", options = list(container = "body")),
                  tipify(colourInput("ma.insig.color", "Insig colour", value = "black"),
                         "Color of insignificant genes.", "right", options = list(container = "body")),
                  tipify(numericInput("ma.sig.opa", label = "Sig opacity:", value = 1, step = 0.05, min = 0),
                         "Opacity of DE genes.", "right", options = list(container = "body"))
                ),
                column(width = 6,
                  numericInput("ma.y", label = "y-axis limits:", value = 5, step = 0.1, min = 0.1),
                  tipify(numericInput("ma.lab.size", label = "Label size:", value = 10, step = 0.5, min = 1),
                         "Font size of gene labels.", "right", options = list(container = "body")),
                  tipify(numericInput("ma.insig.opa", label = "Insig opacity:", value = 1, step = 0.05, min = 0),
                         "Opacity of non-DE genes.", "right", options = list(container = "body")),
                  tipify(numericInput("ma.sig.size", label = "Sig pt size:", value = 5, step = 0.1, min = 1),
                         "Point size of DE genes.", "right", options = list(container = "body")),
                  tipify(numericInput("ma.insig.size", label = "Insig pt size:", value = 3, step = 0.1, min = 0),
                         "Point size of non-DE genes.", "right", options = list(container = "body"))
                )
              ),
              tipify(prettyCheckbox("ma.fcline", label = "Show MAplot FC Threshold", value = TRUE,
                            animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                     "Draw lines for FC threshold.", "right", options = list(container = "body")),
              tipify(prettyCheckbox("ma.hl.counts", label = "Show highlight counts", FALSE, bigger = TRUE,
                             animation = "smooth", status = "success",
                             icon = icon("check"), width = "100%"),
                     "Show count of highlighted genes and genesets.", "right", options = list(container = "body")),
              splitLayout(
                tipify(prettyCheckbox("ma.webgl", label = "Use webGL", TRUE, bigger = TRUE,
                               animation = "smooth", status = "success",
                               icon = icon("check"), width = "100%"),
                       "Plot with webGL. Faster, but sometimes has visual artifacts.", "right", options = list(container = "body")),
                tipify(prettyCheckbox("ma.counts", label = "Show counts", TRUE, bigger = TRUE,
                               animation = "smooth", status = "success",
                               icon = icon("check"), width = "100%"),
                       "Show count of DE and total genes.", "right", options = list(container = "body"))
              ),
              splitLayout(
                tipify(numericInput("ma.webgl.ratio", label = "webGL pixel ratio:", value = 7, step = 0.1, min = 1),
                       "Controls rasterization resolution when webGL is used. Higher is greater resolution, recommend leaving at default.",
                       "right", options = list(container = "body")),
                tipify(numericInput("ma.counts.size", label = "Counts size:", value = 8, step = 0.1, min = 0),
                       "Font size of gene counts.", "right", options = list(container = "body"))
              ),
              circle = FALSE, label = strong("MA-Plot"), status = "danger", size = "lg", icon = icon("gear"),
              width = "300px", tooltip = tooltipOptions(title = "Click to change plot settings")
            ),
            withSpinner(
              jqui_resizable(
                plotlyOutput("ma_plot", height = "500px", width = "550px")
              )
            ),
            div(downloadButton("download_plotly_ma", "Download Interactive MAplot"), align = "left", style="margin-top: 10px;")
          ),
          column(width = 6,
            br(),
            dropdownButton(
              tags$h3("Plot Settings"),
              fluidRow(
                column(width = 6,
                  tipify(colourInput("vol.down.color", "Down-reg colour", value = "#0026ff"),
                         "Color of down-regulated genes.", "right", options = list(container = "body")),
                  tipify(colourInput("vol.up.color", "Up-reg colour", value = "red"),
                         "Color of up-regulated genes.", "right", options = list(container = "body")),
                  tipify(colourInput("vol.insig.color", "Insig colour", value = "#A6A6A6"),
                         "Color of insignificant genes.", "right", options = list(container = "body")),
                  tipify(numericInput("vol.sig.opa", label = "Sig opacity:", value = 1, step = 0.05, min = 0),
                         "Opacity of DE genes.", "right", options = list(container = "body")),
                  tipify(numericInput("vol.sig.size", label = "Sig pt size:", value = 5, step = 0.1, min = 0),
                         "Point size of DE genes.", "right", options = list(container = "body"))
                ),
                column(width = 6,
                  numericInput("vol.x", label = "x-axis limits:", value = 5, step = 0.1, min = 0.1),
                  numericInput("vol.y", label = "y-axis limits:",
                               value = max(-log10(res[[sig.term]][!is.na(res[[sig.term]])])) + 0.1,
                               step = 0.5, min = 1),
                  tipify(numericInput("vol.lab.size", label = "Label size:", value = 10, step = 0.5, min = 1),
                         "Font size of gene labels.", "right", options = list(container = "body")),
                  tipify(numericInput("vol.insig.opa", label = "Insig opacity:", value = 1, step = 0.05, min = 0),
                         "Opacity of non-DE genes.", "right", options = list(container = "body")),
                  tipify(numericInput("vol.insig.size", label = "Insig pt size:", value = 3, step = 0.1, min = 0),
                         "Point size of non-DE genes.", "right", options = list(container = "body")),
                )
              ),
              prettyCheckbox("vol.fcline", label = "Show FC Threshold", value = TRUE,
                             animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
              tipify(prettyCheckbox("vol.sigline", label = "Show Signficance Threshold", value = TRUE,
                             animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                     "Draw line at significance threshold.", "right", options = list(container = "body")),
              tipify(prettyCheckbox("vol.hl.counts", label = "Show Highlight Counts", value = FALSE,
                             animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                     "Show count of highlighted genes and genesets.", "right", options = list(container = "body")),
              splitLayout(
                tipify(prettyCheckbox("vol.counts", label = "Show counts", TRUE, bigger = TRUE,
                               animation = "smooth", status = "success",
                               icon = icon("check"), width = "100%"),
                       "Show count of DE and total genes.", "right", options = list(container = "body")),
                tipify(prettyCheckbox("vol.webgl", label = "Use webGL", TRUE, bigger = TRUE,
                               animation = "smooth", status = "success",
                               icon = icon("check"), width = "100%"),
                       "Plot with webGL. Faster, but sometimes has visual artifacts.", "right", options = list(container = "body"))
              ),
              splitLayout(
                tipify(numericInput("vol.counts.size", label = "Counts size:", value = 8, step = 0.1, min = 0),
                       "Font size of gene counts.", "right", options = list(container = "body")),
                tipify(numericInput("vol.webgl.ratio", label = "webGL pixel ratio:", value = 7, step = 0.1, min = 1),
                       "Controls rasterization resolution when webGL is used. Higher is greater resolution, recommend leaving at default.",
                       "right", options = list(container = "body"))
              ),
              circle = FALSE, label = strong("Volcano Plot"), status = "danger", size = "lg", icon = icon("gear"),
              width = "300px", tooltip = tooltipOptions(title = "Click to change plot settings")
            ),
            withSpinner(jqui_resizable(plotlyOutput("volcano_plot", height = "500px", width = "550px"))),
            div(downloadButton("download_plotly_volc", "Download Interactive Volcano plot"), align = "left", style="margin-top: 10px;")
          )
        )
      ),
      tabPanel("Full DE Results Table",
        div(DT::dataTableOutput("res_table_full"), style = "font-size:80%; margin:10px;")
      ),
      tabPanel("Sample Metadata (Filtering)",
        div(DT::dataTableOutput("metadata"), style = "font-size:80%; margin:10px;")
      )
    )
  )

  # Side bar contains settings for certain cutoffs to select significant genes.
  ui <- dashboardPage(
    dashboardHeader(disable = TRUE),
    dashboardSidebar(disable = TRUE),
    dashboardBody(
      tags$head(
        # Note the wrapping of the string in HTML()
        tags$style(HTML("
          .panel-body {
            padding: 5px;
          }
          .form-group {
            margin-bottom: 5px;
          }
          .well {
            padding: 5px;
            margin-bottom: 10px;
          }
          label {
            font-size: 80%;
          }
          .form-control, .selectize-input{
            padding-bottom: 2px !important;
            padding-top: 2px !important;
            font-size: 10px;
            height: 24px;
          }
        "))
      ),
      useShinyjs(),
      shinyDashboardThemes(
        theme = "onenote"
      ),
      sidebarLayout(
        sidebarPanel(
          width = 2,
          tags$label(HTML(qq("Comparison: <code style='font-weight:normal; font-size: 10px;'>@{paste(coef, collapse = ' ')}</code>")),
                     class = "shiny-input-container", style = "font-size:1.2em;"),
          hidden(div(id = "mres", tipify(selectInput("res.select", NULL, choices = names(res.list)),
                 "Results to view from provided set.", "right", options = list(container = "body")))),
          hr(style="margin:2px; background-color: #737373;"),
          bsCollapse(open = "settings",
            bsCollapsePanel(title = span(icon("plus"), "Plot Settings"), value = "settings", style = "info",
              tipify(numericInput("sig.thresh", label = "Significance threshold:", value = 0.05, step = 0.001, min = 0.0001),
                     "Significance threshold to consider a gene DE.", "right", options = list(container = "body")),
              tipify(selectInput("sig.term", label = "Significance term:", choices = colnames(res),
                          selected = ifelse("padj" %in% colnames(res), "padj", "svalue")),
                     "Significance term to use for DE filtering.", "right", options = list(container = "body")),
              tipify(numericInput("base_mean", label = "Minimal baseMean:", value = 0, step = 1),
                     "Minimal baseMean (expression) required to consider a gene DE.", "right", options = list(container = "body")),
              tipify(numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 0, step = 0.1, min = 0),
                     "log2 fold change magnitude threshold required to consider a gene DE.", "right", options = list(container = "body")),
              tipify(numericInput("row.km", label = "Row k-means groups:", value = 2, step = 1),
                     "Number of groups to break heatmap into via  k-means clustering on rows.", "right", options = list(container = "body")),
              tipify(numericInput("col.km", label = "Column k-means groups:", value = 0, step = 1),
                     "Number of groups to break heatmap into via k-means clustering on columns.", "right", options = list(container = "body")),
              tipify(pickerInput("anno.vars", "Annotate Heatmap by:", choices = c("", names(SummarizedExperiment::colData(dds))),
                          multiple = TRUE, options = list(`live-search` = TRUE, `actions-box` = TRUE),
                          selected = annot.by),
                     "Sample metadata columns used for column annotations.", "right", options = list(container = "body"))
            ),
            bsCollapsePanel(title = span(icon("plus"), "Highlight Gene(sets)"), value = "genesets", style = "info",
              textAreaInput("hl.genes", "Highlight Genes:", value = "", rows = 4,
                            placeholder = "Enter space, comma, or newline delimited genes"),
              pickerInput("hl.genesets", "Highlight Genesets:", choices = c("", names(genesets)),
                          multiple = TRUE, options = list(`live-search` = TRUE,
                                                          `actions-box` = TRUE)),
              fluidRow(
                column(6,
                       tipify(numericInput("hl.genes.opa", label = "Genes opacity:", value = 1, step = 0.05, min = 0),
                              "Opacity of highlighted genes.", "right", options = list(container = "body")),
                       tipify(numericInput("hl.genes.size", label = "Genes pt size:", value = 7, step = 0.1, min = 0),
                              "Point size of highlighted genes.", "right", options = list(container = "body")),
                       tipify(numericInput("hl.genes.lw", label = "Genes border width:", value = 0.5, step = 0.05, min = 0),
                              "Width of border for highlighted genes.", "right", options = list(container = "body")),
                       tipify(colourInput("hl.genes.col", "Genes color:", value = "#FFFF19"),
                              "Color of genes to highlight.", "right", options = list(container = "body")),
                       tipify(colourInput("hl.genes.lcol", "Genes border:", value = "#000000"),
                              "Border color of genes to highlight.", "right", options = list(container = "body"))
                ),
                column(6,
                       tipify(numericInput("hl.genesets.opa", label = "Sets opacity:", value = 1, step = 0.05, min = 0),
                              "Opacity of highlighted genesets.", "right", options = list(container = "body")),
                       tipify(numericInput("hl.genesets.size", label = "Sets pt size:", value = 7, step = 0.1, min = 0),
                              "Point size of highlighted genesets.", "right", options = list(container = "body")),
                       tipify(numericInput("hl.genesets.lw", label = "Sets border width:", value = 0.5, step = 0.05, min = 0),
                              "Width of border for highlighted genesets.", "right", options = list(container = "body")),
                       tipify(colourInput("hl.genesets.col", "Sets color:", value = "#38FFF2"),
                              "Color of genesets to highlight.", "right", options = list(container = "body")),
                       tipify(colourInput("hl.genesets.lcol", "Sets border:", value = "#000000"),
                              "Border color of genesets to highlight.", "right", options = list(container = "body")))
              )
            )
          ),
          div(actionButton("update", label = "Update Plots"), align = "center")
        ),
        body
      )
    )
  )

  server <- function(input, output, session) {

    if (multi.res) {
      shinyjs::show("mres")
    }

    # Get annotations. If none provided, use design variables.
    anno <- reactiveVal()

    # Used to hold plots for download.
    plot_store <- reactiveValues()

    observeEvent(c(input$anno.vars,
                   input$res.select,
                   input$metadata_rows_all,
                   input$update), {

      annos <- SummarizedExperiment::colData(dds)

      if (!is.null(input$anno.vars)) {
        annos <- annos[, input$anno.vars, drop = FALSE]
      } else {
        annos <- NULL
      }

      if (!is.null(annos)) {
        l <- sapply(annos, function(x) (is.factor(x) || is.character(x))) | sapply(annos, function(x) length(unique(x)) == 1)
        annos <- annos[, l, drop = FALSE]
      }

      if (ncol(annos) == 0) {
        annos <- NULL
      }

      # Filter samples.
      if (!is.null(input$metadata_rows_all)) {
        annos <- annos[input$metadata_rows_all,]
      }

      anno(annos)
    })

    # Keep track of which genes have been clicked
    genes <- reactiveValues(ma = NULL, volc = NULL)

    ress <- reactiveVal({res})

    # Get the selected results tables.
    observeEvent(input$res.select, {
      req(ress, anno)
      ress(res.list[[input$res.select]])

      if (!is.null(input$metadata_rows_all)) {
        mat <- mat[, input$metadata_rows_all]
      }

      pdf(NULL)
      ht <- .make_heatmap(mat, ress(), anno(), baseMean_col_fun, log2fc_col_fun, sig.term = input$sig.term,
                          sig.thresh = as.numeric(input$sig.thresh), base_mean = input$base_mean, log2fc = input$log2fc,
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
    }, ignoreInit = TRUE)

    # A self-defined action to respond to heatmap selections to show them in a table.
    .brush_action <- function(df, input, output, session) {

      row_index <- unique(unlist(df$row_index))
      selected <- env$row_index[row_index]

      cnames <- c(sig.term)
      if (!is.null(swap.rownames)) {
        cnames <- c(cnames, "ORIGINAL_ROWS")
      }

      output[["res_table"]] <- DT::renderDataTable({

        DT::formatRound(DT::datatable(as.data.frame(ress()[selected, c("baseMean", "log2FoldChange", cnames)]),
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
        cnames <- c(sig.term)
        if (!is.null(swap.rownames)) {
          cnames <- c(cnames, "ORIGINAL_ROWS")
        }

        df <- df[, c("baseMean", "log2FoldChange", "lfcSE", cnames)]
        DT::formatRound(DT::datatable(as.data.frame(ress()[selected, c("baseMean", "log2FoldChange", cnames)]),
                                      rownames = TRUE, options = list(lengthMenu = c(5, 10, 25),
                                                                      pageLength = 20)),
                        columns = 1:3, digits = 5) %>%
          DT::formatStyle(0, target = "row", lineHeight = '40%')
      })
    }

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
      input$update

      plot_store$ma_plot <- .make_maplot(res = ress(),
                   ylim = isolate(input$ma.y),
                   fc.thresh = isolate(input$log2fc),
                   fc.lines = isolate(input$ma.fcline),
                   sig.thresh = isolate(input$sig.thresh),
                   basemean.thresh = isolate(input$base_mean),
                   h.id = h.id,
                   sig.term = isolate(input$sig.term),
                   gs = genes$ma,
                   up.color = isolate(input$ma.up.color),
                   down.color = isolate(input$ma.down.color),
                   insig.color = isolate(input$ma.insig.color),
                   sig.opacity = isolate(input$ma.sig.opa),
                   insig.opacity = isolate(input$ma.insig.opa),
                   sig.size = isolate(input$ma.sig.size),
                   insig.size = isolate(input$ma.insig.size),
                   label.size = isolate(input$ma.lab.size),
                   webgl = isolate(input$ma.webgl),
                   webgl.ratio = isolate(input$ma.webgl.ratio),
                   show.counts = isolate(input$ma.counts),
                   show.hl.counts = isolate(input$ma.hl.counts),
                   counts.size = isolate(input$ma.counts.size),
                   highlight.genesets = isolate(input$hl.genesets),
                   highlight.genes = isolate(input$hl.genes),
                   genesets = genesets,
                   highlight.genes.color = isolate(input$hl.genes.col),
                   highlight.genes.size = isolate(input$hl.genes.size),
                   highlight.genes.opac = isolate(input$hl.genes.opa),
                   highlight.genes.linecolor = isolate(input$hl.genes.lcol),
                   highlight.genes.linewidth = isolate(input$hl.genes.lw),
                   highlight.genesets.color = isolate(input$hl.genesets.col),
                   highlight.genesets.size = isolate(input$hl.genesets.size),
                   highlight.genesets.opac = isolate(input$hl.genesets.opa),
                   highlight.genesets.linecolor = isolate(input$hl.genesets.lcol),
                   highlight.genesets.linewidth = isolate(input$hl.genesets.lw))

      plot_store$ma_plot
    })

    output$volcano_plot <- renderPlotly({
      req(genes)
      input$update

      plot_store$volcano_plot <- .make_volcano(res = ress(),
                    xlim = isolate(input$vol.x),
                    ylim = isolate(input$vol.y),
                    fc.thresh = isolate(input$log2fc),
                    fc.lines = isolate(input$vol.fcline),
                    sig.thresh = isolate(input$sig.thresh),
                    sig.line = isolate(input$vol.sigline),
                    basemean.thresh = isolate(input$base_mean),
                    h.id = h.id,
                    sig.term = isolate(input$sig.term),
                    lfc.term = "log2FoldChange",
                    feat.term = "rows",
                    fs = genes$volc,
                    up.color = isolate(input$vol.up.color),
                    down.color = isolate(input$vol.down.color),
                    insig.color = isolate(input$vol.insig.color),
                    sig.opacity = isolate(input$vol.sig.opa),
                    insig.opacity = isolate(input$vol.insig.opa),
                    sig.size = isolate(input$vol.sig.size),
                    insig.size = isolate(input$vol.insig.size),
                    label.size = isolate(input$vol.lab.size),
                    webgl = isolate(input$vol.webgl),
                    webgl.ratio = isolate(input$vol.webgl.ratio),
                    show.counts = isolate(input$vol.counts),
                    show.hl.counts = isolate(input$vol.hl.counts),
                    counts.size = isolate(input$vol.counts.size),
                    highlight.featsets = isolate(input$hl.genesets),
                    highlight.feats = isolate(input$hl.genes),
                    featsets = genesets,
                    highlight.feats.color = isolate(input$hl.genes.col),
                    highlight.feats.size = isolate(input$hl.genes.size),
                    highlight.feats.opac = isolate(input$hl.genes.opa),
                    highlight.feats.linecolor = isolate(input$hl.genes.lcol),
                    highlight.feats.linewidth = isolate(input$hl.genes.lw),
                    highlight.featsets.color = isolate(input$hl.genesets.col),
                    highlight.featsets.size = isolate(input$hl.genesets.size),
                    highlight.featsets.opac = isolate(input$hl.genesets.opa),
                    highlight.featsets.linecolor = isolate(input$hl.genesets.lcol),
                    highlight.featsets.linewidth = isolate(input$hl.genesets.lw))

      plot_store$volcano_plot
    })

    output[["res_table_full"]] <- DT::renderDT(server = FALSE, {
      req(ress)
      df <- as.data.frame(ress())

      # Collect proper output columns.
      cnames <- NULL
      if (sig.term == "svalue") {
        cnames <- "svalue"
      } else if ("stat" %in% colnames(df)) {
        cnames <- c("padj", "pvalue", "stat")
      } else {
        cnames <- c("padj", "pvalue")
      }

      # Can't tack on to cnames because of rounding.
      snames <- NULL
      if (!is.null(swap.rownames)) {
        snames <- "ORIGINAL_ROWS"
      }

      df <- df[, c("baseMean", "log2FoldChange", "lfcSE", cnames, snames)]

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

    # Metadata table.
    output$metadata <- DT::renderDT(server = FALSE, {
      df <- as.data.frame(SummarizedExperiment::colData(dds))
      DT::datatable(df,
                    filter = "top",
                    extensions = c("Buttons", "Scroller"),
                    options = list(
                      search = list(regex = TRUE),
                      lengthMenu = list(c(10, 25, 50, -1), c("10", "25", "50", "all")),
                      dom = 'Blfrtip',
                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                      scrollX = TRUE,
                      deferRender = TRUE,
                      scrollY = 400,
                      scroller = TRUE)
      )
    })

    # Only remake heatmap on button click.
    observeEvent(input$update, {
      req(ress, anno)

      if (is.null(input$anno.vars)) {
        anno(NULL)
      }

      if (!is.null(input$metadata_rows_all)) {
        mat <- mat[,input$metadata_rows_all]
      }

      pdf(NULL)
      ht <- .make_heatmap(mat, ress(), anno(), baseMean_col_fun, log2fc_col_fun, sig.term = input$sig.term,
                          sig.thresh = as.numeric(input$sig.thresh), base_mean = input$base_mean, log2fc = input$log2fc,
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

    # Download interactive plots as html.
    output$download_plotly_volc <- downloadHandler(
      filename = function() {
        paste("volcanoplot-", Sys.Date(), ".html", sep = "")
      },
      content = function(file) {
        # export plotly html widget as a temp file to download.
        saveWidget(jqui_resizable(plot_store$volcano_plot),
                   file, selfcontained = TRUE)
      }
    )

    output$download_plotly_ma <- downloadHandler(
      filename = function() {
        paste("maplot-", Sys.Date(), ".html", sep = "")
      },
      content = function(file) {
        # export plotly html widget as a temp file to download.
        saveWidget(jqui_resizable(plot_store$ma_plot),
                   file, selfcontained = TRUE)
      }
    )
  }

  shinyApp(ui, server, options = list(height = height))
}
