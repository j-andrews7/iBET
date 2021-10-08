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
#' @importFrom shinycustomloader withLoader
#' @importFrom shinyjqui jqui_resizable
#' @importFrom colourpicker colourInput
#' @importFrom utils combn
#' @importFrom stats lm fitted.values
#'
#' @param res A named list of data.frames containing differential expression analysis results.
#' @param sig.col String for the column name of the significance value, e.g. "padj". If not provided,
#'   the function will search for commonly used values ("padj", "FDR", "svalue", "adj.P.Val") in the column names.
#' @param sig.thresh Number to be used as the significance threshold. Adjustable within the app.
#' @param lfc.col String for the column name of the log2 fold change column, e.g. "log2FC". If not provided,
#'   the function will search for commonly used values c("log2FoldChange", "logFC") in the column names.
#' @param gene.col String for the column name containing the gene identifier. If not provided, rownames will
#'   be used.
#' @param expr.col Optional string for the column name containing average expression. If not provided,
#'   the function will search for commonly used values ("baseMean", "logCPM", "AveExpr") in the column names.
#' @param height Number indicating height of app in pixels.
#'
#' @return A Shiny app containing scatter plots comparing all combinations of the DE results to each other,
#' along with a line of best fit and correlation testing.
#'
#'
#' @author Jared Andrews
#' @export
shinyDECorr <- function(res, sig.col = NULL, sig.thresh = 0.05, lfc.col = NULL,
                        gene.col = NULL, expr.col = NULL, height = 800) {

  # Parameter validation.
  if (is.null(names(res))) {
    stop("Results list must be named")
  }

  if (length(res) < 2) {
    stop("Results list must contain at least 2 elements")
  } else if (length(res) > 4) {
    stop("A maximum of 4 DE results can be provided")
  }

  res.names <- colnames(res[[1]])

  if (is.null(sig.col)) {
    if (!any(res.names %in% c("padj", "FDR", "svalue", "adj.P.Val"))) {
      stop("Cannot determine significance column, please provide the column name to sig.col")
    } else {
      sig.col <- res.names[res.names %in% c("padj", "FDR", "svalue", "adj.P.Val")]
    }
  }

  if (is.null(lfc.col)) {
    if (!any(res.names %in% c("log2FoldChange", "logFC"))) {
      stop("Cannot determine fold change column, please provide the column name to lfc.col")
    } else {
      lfc.col <- res.names[res.names %in% c("log2FoldChange", "logFC")]
    }
  }

  if (is.null(expr.col)) {
    if (!any(res.names %in% c("baseMean", "logCPM", "AveExpr"))) {
      message("Cannot determine average expression column")
    } else {
      expr.col <- res.names[res.names %in% c("baseMean", "logCPM", "AveExpr")]
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
      shinyDashboardThemes(theme = "onenote"),
      sidebarLayout(
        sidebarPanel(
          width = 2,
          h3("Plot Settings"),
          numericInput("sig", label = "Significance threshold:", value = 0.05, step = 0.001, min = 0.0001),
          numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 0, step = 0.1, min = 0),
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
          hr(),
          h4("Point Color"),
          splitLayout(
            colourInput("comp1.sig", "X-axis Signif", value = "#E69F00"),
            colourInput("comp2.sig", "Y-axis Signif", value = "#56B4E9")
          ),
          splitLayout(
            colourInput("both.sig", "Both Signif", value = "#009E73"),
            colourInput("insig.color", "Insignificant", value = "#666666")
          )
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
          withLoader(
            jqui_resizable(
              plotlyOutput(row1[x], height = "350px", width = "350px")
            ),
            type = "html", loader = "dnaspin"
          )
        )
      })

      # ecessary for the list of
      # items to display properly
      do.call(tagList, row1_plots)
    })

    if (!is.null(row2)) {

      output$row2 <- renderUI({

        row2_plots <- lapply(1:length(row2), function(x) {
          column(width = 4,
            withLoader(
              jqui_resizable(
                plotlyOutput(row2[x], height = "350px", width = "350px")
              ),
              type = "html", loader = "dnaspin"
            )
          )
        })

        # Necessary for the list of items to display properly
        do.call(tagList, row2_plots)
      })
    } else {
      output$row2 <- renderUI({div()})
    }

    for (n in colnames(res.comb)) {
      local({
        my_n <- n
        df1 <- res.comb[[1, my_n]]
        df2 <- res.comb[[2, my_n]]

        # Make plots. Has to be plot_ly made to add the labels. Regression line still confusing.
        output[[my_n]] <- renderPlotly({
          req(genes)
          .make_xyplot(res[df1], res[df2], sig.col = sig.col, lfc.col = lfc.col,
                       sig.thresh = input$sig,
                       lfc.thresh = input$log2fc, gene.col = gene.col, source = my_n,
                       expr.col = expr.col, regr = input$draw.reg, genes.labeled = genes[[my_n]],
                       res1.color = input$comp1.sig, res2.color = input$comp2.sig,
                       both.color = input$both.sig, insig.color = input$insig.color,
                       xlim = input$xlim, ylim = input$ylim, show = input$show)
        })
      })
    }
  }

  shinyApp(ui, server, options = list(height = height))
}
