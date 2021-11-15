#' Create an interactive Shiny app for correlation of differential expression results
#'
#' This shiny app generates scatter plots for every combination of DE results fed to it. It is useful for
#' comparing relative differences in differential expression for common conditions between different
#' backgrounds, groups, or settings.
#'
#' Gene labels can be added to the plots by clicking a point. The labels can also be dragged around,
#' though adding labels to a plot will reset the positions, so it's recommended to add all labels prior to re-positioning them.
#'
#' @details
#' Comparisons will be limited to shared genes.
#'
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom plotly ggplotly plotlyOutput renderPlotly toWebGL plot_ly layout add_annotations config toRGB event_data
#' @import ggplot2
#' @import shinydashboard
#' @import dashboardthemes
#' @importFrom shinyWidgets prettyCheckbox pickerInput tooltipOptions
#' @importFrom shinyjqui jqui_resizable
#' @importFrom colourpicker colourInput
#' @importFrom utils combn
#' @importFrom stats lm fitted.values
#' @importFrom shinyBS bsCollapse bsCollapsePanel
#'
#' @param res A named list of data.frames containing differential expression analysis results.
#' @param sig.col String for the column name of the significance value, e.g. "padj". If not provided,
#'   the function will search for commonly used values ("padj", "FDR", "svalue", "adj.P.Val") in the column names.
#' @param sig.thresh Number to be used as the significance threshold. Adjustable within the app.
#' @param lfc.col String for the column name of the log2 fold change column, e.g. "log2FC". If not provided,
#'   the function will search for commonly used values ("log2FoldChange", "logFC", "LFC") in the column names.
#' @param gene.col String for the column name containing the gene identifier. If not provided, rownames will
#'   be used.
#' @param expr.col Optional string for the column name containing average expression. If not provided,
#'   the function will search for commonly used values ("baseMean", "logCPM", "AveExpr") in the column names.
#' @param genesets Optional named list containing genesets that can be interactively highlighted on the plots.
#'   The elements of the list should each be a geneset with gene identifiers matching those used in the results.
#' @param height Number indicating height of app in pixels.
#'
#' @return A Shiny app containing scatter plots comparing all combinations of the DE results to each other,
#' along with a line of best fit, correlation testing, gene and geneset highlighting, and movable annotations.
#'
#'
#' @author Jared Andrews
#' @export
shinyDECorr <- function(res, sig.col = NULL, sig.thresh = 0.05, lfc.col = NULL,
                        gene.col = NULL, expr.col = NULL, genesets = NULL, height = 800) {

  # Parameter validation.
  if (is.null(names(res))) {
    stop("Results list must be named")
  }

  if (length(res) < 2) {
    stop("Results list must contain at least 2 elements")
  } else if (length(res) > 4) {
    stop("A maximum of 4 DE results can be provided")
  }

  if (!is.null(genesets)) {
    if (is.null(names(genesets))) {
      stop("Genesets list must be named")
    } else if (!is(genesets, "list")) {
      stop("Genesets must be provided as a named list")
    }
  }

  # Get all combinations and rows needed.
  res.comb <- combn(names(res), 2)
  colnames(res.comb) <- apply(res.comb, 2, paste0, collapse = "")

  if (ncol(res.comb) > 3) {
    row1 <- colnames(res.comb)[1:3]
    row2 <- colnames(res.comb)[4:ncol(res.comb)]
  } else {
    row1 <- colnames(res.comb)[1:ncol(res.comb)]
    row2 <- NULL
  }

  body <- mainPanel(width = 10,
    fluidRow(
      uiOutput("row1")
    ),
    fluidRow(
      uiOutput("row2")
    )
  )

  # For label tracking.
  genes.init <- list()
  for (n in colnames(res.comb)) {
    genes.init[[n]] <- NULL
  }

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
        "))
      ),
      shinyDashboardThemes(theme = "onenote"),
      sidebarLayout(
        sidebarPanel(
          width = 2,
          bsCollapse(open = "settings",
            bsCollapsePanel(title = span(icon("plus"), "Plot Settings"), value = "settings", style = "info",
              splitLayout(
                numericInput("sig", label = "Sig. threshold:", value = 0.05, step = 0.001, min = 0.0001),
                numericInput("log2fc", label = "log2FC theshold:", value = 0, step = 0.1, min = 0)
              ),
              splitLayout(
                numericInput("ylim", label = "Y-axis limit:", value = 10, step = 0.1, min = 0),
                numericInput("xlim", label = "X-axis limit:", value = 10, step = 0.1, min = 0)
              ),
              pickerInput("show", label = "Display genes:", choices = c("Both Significant", "X-axis Significant",
                                                                       "Y-axis Significant", "Not Significant"),
                          selected = c("Both Significant", "X-axis Significant",
                                       "Y-axis Significant"), multiple = TRUE),
              prettyCheckbox("draw.reg", strong("Draw regression line"), TRUE, bigger = TRUE,
                             animation = "smooth", status = "success",
                             icon = icon("check"), width = "100%"),
              prettyCheckbox("webgl", strong("Use webGL"), TRUE, bigger = TRUE,
                             animation = "smooth", status = "success",
                             icon = icon("check"), width = "100%"),
              numericInput("webgl.ratio", label = "webGL pixel ratio:", value = 7, step = 0.1, min = 1),
              prettyCheckbox("counts", strong("Show counts"), TRUE, bigger = TRUE,
                             animation = "smooth", status = "success",
                             icon = icon("check"), width = "100%"),
              prettyCheckbox("hl.counts", strong("Show highlight counts"), FALSE, bigger = TRUE,
                             animation = "smooth", status = "success",
                             icon = icon("check"), width = "100%"),
              splitLayout(
                numericInput("aggr.size", label = "Corr stats size:", value = 8, step = 0.1, min = 0),
                numericInput("counts.size", label = "Counts size:", value = 8, step = 0.1, min = 0)
              )
            ),
            bsCollapsePanel(title = span(icon("plus"), "Point Aesthetics"), style = "info",
              numericInput("lab.size", label = "Label Size:", value = 10, step = 0.5, min = 1),
              fluidRow(
                column(6,
                       numericInput("x.opa", label = "X-sig opacity:", value = 1, step = 0.05, min = 0),
                       numericInput("y.opa", label = "Y-sig opacity:", value = 1, step = 0.05, min = 0),
                       numericInput("x.size", label = "X-sig pt size:", value = 5, step = 0.1, min = 0),
                       numericInput("y.size", label = "Y-sig pt size:", value = 5, step = 0.1, min = 0),
                       colourInput("comp1.sig", "X-axis Signif", value = "#E69F00"),
                       colourInput("both.sig", "Both Signif", value = "#009E73")),
                column(6,
                       numericInput("both.opa", label = "Both opacity:", value = 1, step = 0.05, min = 0),
                       numericInput("insig.opa", label = "Insig opacity:", value = 1, step = 0.05, min = 0),
                       numericInput("both.size", label = "Both pt size:", value = 5, step = 0.1, min = 0),
                       numericInput("insig.size", label = "Insig pt size:", value = 3, step = 0.1, min = 0),
                       colourInput("comp2.sig", "Y-axis Signif", value = "#BC57EB"),
                       colourInput("insig.color", "Insignificant", value = "#666666"))
              )
            ),
            bsCollapsePanel(title = span(icon("plus"), "Highlight Gene(sets)"), style = "info",
              textAreaInput("hl.genes", "Highlight Genes:", value = "", rows = 3,
                            placeholder = "Enter space, comma, or newline delimited genes"),
              pickerInput("hl.genesets", "Highlight Genesets:", choices = c("", names(genesets)),
                          multiple = TRUE, options = list(`live-search` = TRUE,
                                                          `actions-box` = TRUE)),
              fluidRow(
                column(6,
                  numericInput("hl.genes.opa", label = "Genes opacity:", value = 1, step = 0.05, min = 0),
                  numericInput("hl.genes.size", label = "Genes pt size:", value = 7, step = 0.1, min = 0),
                  numericInput("hl.genes.lw", label = "Genes border width:", value = 0.5, step = 0.05, min = 0),
                  colourInput("hl.genes.col", "Genes color:", value = "#FFFF19"),
                  colourInput("hl.genes.lcol", "Genes border:", value = "#000000")),
                column(6,
                  numericInput("hl.genesets.opa", label = "Sets opacity:", value = 1, step = 0.05, min = 0),
                  numericInput("hl.genesets.size", label = "Sets pt size:", value = 7, step = 0.1, min = 0),
                  numericInput("hl.genesets.lw", label = "Sets border width:", value = 0.5, step = 0.05, min = 0),
                  colourInput("hl.genesets.col", "Sets color:", value = "#38FFF2"),
                  colourInput("hl.genesets.lcol", "Sets border:", value = "#000000"))
              )
            )
          ),
          div(actionButton("update", "Update Plots"), align = "center")
        ),
        body
      )
    )
  )

  server <- function(input, output, session) {

    # Keep track of which genes have been clicked
    genes <- do.call("reactiveValues", genes.init)

    # On click, the key field of the event data contains the gene symbol
    # Add that gene to the set of all "selected" genes
    lapply(1:length(colnames(res.comb)), FUN = function(x) {
      n <- colnames(res.comb)[x]
      observeEvent(event_data("plotly_click", source = n), {
        gene <- event_data("plotly_click", source = n)
        gene_old_new <- rbind(genes[[n]], gene)
        keep <- gene_old_new[gene_old_new$customdata %in% names(which(table(gene_old_new$customdata)==1)),]

        if (nrow(keep) == 0) {
          genes[[n]] <- NULL
        } else {
          genes[[n]] <- keep
        }
      })
    })

    # clear the set of genes when a double-click occurs
    lapply(1:length(colnames(res.comb)), FUN = function(x) {
      n <- colnames(res.comb)[x]
      observeEvent(event_data("plotly_doubleclick", source = n), {
        genes[[n]] <- NULL
      })
    })

    output$row1 <- renderUI({
      req(genes)
      # dynamically allocate rows/columns based on number of plots
      row1_plots <- lapply(1:length(row1), function(x) {
        column(width = 4,
          jqui_resizable(
            plotlyOutput(row1[x], height = "350px", width = "350px")
          )
        )
      })

      # Necessary for the list of items to display properly
      do.call(tagList, row1_plots)
    })

    if (!is.null(row2)) {

      output$row2 <- renderUI({

        row2_plots <- lapply(1:length(row2), function(x) {
          column(width = 4,
            jqui_resizable(
              plotlyOutput(row2[x], height = "350px", width = "350px")
            )
          )
        })

        # Necessary for the list of items to display properly
        do.call(tagList, row2_plots)
      })
    } else {
      output$row2 <- renderUI({div()})
    }

    # Iteratively make plots.
    for (n in colnames(res.comb)) {
      local({
        my_n <- n
        df1 <- res.comb[[1, my_n]]
        df2 <- res.comb[[2, my_n]]

        df.vars <- .get_plot_vars(res[df1], res[df2], sig.col = sig.col,
                                  lfc.col = lfc.col, expr.col = expr.col)

        output[[my_n]] <- renderPlotly({
          req(genes)
          input$update
          .make_xyplot(res[df1], res[df2],
                       df.vars = df.vars,
                       sig.thresh = isolate(input$sig),
                       lfc.thresh = isolate(input$log2fc),
                       gene.col = gene.col,
                       source = my_n,
                       regr = isolate(input$draw.reg),
                       genes.labeled = genes[[my_n]],
                       res1.color = isolate(input$comp1.sig),
                       res2.color = isolate(input$comp2.sig),
                       both.color = isolate(input$both.sig),
                       insig.color = isolate(input$insig.color),
                       xlim = isolate(input$xlim),
                       ylim = isolate(input$ylim),
                       show = isolate(input$show),
                       label.size = isolate(input$lab.size),
                       webgl = isolate(input$webgl),
                       webgl.ratio = isolate(input$webgl.ratio),
                       show.counts = isolate(input$counts),
                       counts.size = isolate(input$counts.size),
                       show.hl.counts = isolate(input$hl.counts),
                       aggr.size = isolate(input$aggr.size),
                       res1.size = isolate(input$x.size),
                       res2.size = isolate(input$y.size),
                       both.size = isolate(input$both.size),
                       insig.size = isolate(input$insig.size),
                       res1.opac = isolate(input$x.opa),
                       res2.opac = isolate(input$y.opa),
                       both.opac = isolate(input$both.opa),
                       insig.opac = isolate(input$insig.opa),
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
        })
      })
    }
  }

  shinyApp(ui, server, options = list(height = height))
}
