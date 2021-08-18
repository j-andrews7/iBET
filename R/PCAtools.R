#' Create an interactive Shiny app for PCAtools
#'
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @import DT
#' @import PCAtools
#' @importFrom plotly ggplotly plotlyOutput renderPlotly toWebGL
#' @import ggplot2
#' @import shinydashboard
#' @import dashboardthemes
#' @importFrom dittoSeq dittoColors
#' @importFrom shinyWidgets prettyCheckbox
#' @importFrom shinycustomloader withLoader
#' @importFrom shinyjqui jqui_resizable
#' @importFrom matrixStats rowVars
#'
#' @param mat A matrix with features as rows and samples as columns.
#' @param metadata A dataframe containing sample metadata. The rownames must match the column names of the matrix.
#' @param annot.by A string or character vector containing the names of sample metadata variables to be used for hover text.
#' @param color.by A string containing the name of a sample metadata variable to be used to color points.
#' @param shape.by A string containing the name of a sample metadata variable to be used to shape points.
#' @param height Number indicating height of app in pixels.
#' @inheritParams PCAtools::pca
#'
#' @return A Shiny app wrapped around most PCAtools functions and plots.
#'
#'
#' @seealso
#' \code{\link[PCAtools]{pca}}, \code{\link[PCAtools]{screeplot}}, \code{\link[PCAtools]{biplot}}, \code{\link[PCAtools]{pairsplot}}.
#'
#' @author Jared Andrews
#' @export
shinyPCAtools <- function(mat, metadata, removeVar = 0.3, scale = TRUE,
                          center = TRUE, color.by = NULL, shape.by = NULL,
                          annot.by = NULL, height = 800) {

  ui <- dashboardPage(
    dashboardHeader(disable = TRUE),
    dashboardSidebar(disable = TRUE),
    dashboardBody(
      shinyDashboardThemes(
        theme = "onenote"
      ),
      sidebarLayout(
        sidebarPanel(width = 2,
                     h3("Plot Parameters"),
                     uiOutput("pca.comps"),
                     selectInput("color", "Color by:", choices = c("", colnames(metadata)), selected = color.by),
                     selectInput("shape", "Shape by:", choices = c("", colnames(metadata)), selected = shape.by),
                     prettyCheckbox("meta.filt", strong("Limit via metadata table"), TRUE, bigger = TRUE,
                                    animation = "smooth", status = "success",
                                    icon = icon("check"), width = "100%"),
                     h3("PCA Parameters"),
                     numericInput("var.remove", "Remove this proportion of genes based on low variance:",
                                  min = 0, max = 1, step = 0.01, value = removeVar),
                     prettyCheckbox("center", strong("Center data"), center, bigger = TRUE,
                                    animation = "smooth", status = "success",
                                    icon = icon("check"), width = "100%"),
                     prettyCheckbox("scale", strong("Scale data"), scale, bigger = TRUE,
                                    animation = "smooth", status = "success",
                                    icon = icon("check"), width = "100%"),
                     prettyCheckbox("twod", strong("Limit to 2D"), FALSE, bigger = TRUE,
                                    animation = "smooth", status = "success",
                                    icon = icon("check"), width = "100%"),
        ),
        mainPanel(width = 10,
                  tabsetPanel(
                    tabPanel("Plot", div(plotlyOutput("main.plot"), align = "center", style = "height:700px;")),
                    tabPanel("Metadata", div(DTOutput("metadata"), style = "font-size:80%"))
                  )
        )
      )
    )
  )

  server = function(input, output, session) {

    pc <- reactive({
      matty <- mat
      meta <- metadata

      if (!is.null(input$metadata_rows_all) & input$meta.filt) {
        matty <- matty[,input$metadata_rows_all]
        meta <- metadata[input$metadata_rows_all,]
      }

      # Remove features with no variance.
      matty <- matty[rowVars(matty) > 0,]

      pca(matty, metadata = meta, removeVar = input$var.remove, scale = input$scale, center = input$center)
    })

    # Populate UI with all PCs.
    output$pca.comps <- renderUI({
      req(pc)
      pcs <- pc()

      tagList(
        selectInput("dim1", "Dim1:", choices = pcs$components, selected = "PC1"),
        selectInput("dim2", "Dim2:", choices = pcs$components, selected = "PC2"),
        selectInput("dim3", "Dim3:", choices = pcs$components, selected = "PC3")
      )
    })

    output$metadata <- renderDT({
      DT::datatable(as.data.frame(metadata),
                    rownames = FALSE,
                    filter = "top",
                    extensions = c("Buttons", "Scroller"),
                    options = list(
                      search = list(regex = TRUE),
                      lengthMenu = list(c(10, 25, 50, -1), c("10", "25", "50", "all")),
                      dom = 'Blfrtip',
                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                      scrollX = TRUE,
                      deferRender = TRUE,
                      scrollY = 600,
                      scroller = TRUE)
      ) %>% DT::formatStyle(0, target = "row", lineHeight = '40%')
    })

    output$main.plot = renderPlotly({
      req(pc)

      pc.res <- pc()

      pl.cols <- NULL
      pl.shapes <- NULL
      pl.col <- "black"
      hov.text <- NULL

      # Get marker aesthetics mappings.
      # Drop unused factor levels if possible.
      if (input$color != "") {
        pl.cols <- pc.res$metadata[,input$color, drop = TRUE]
        if (is.factor(pl.cols)) {
          pl.cols <- droplevels(pl.cols)
        }
        pl.col <- dittoColors()[seq_along(unique(pc.res$metadata[,input$color, drop = TRUE]))]
      }

      if (input$shape != "") {
        pl.shapes <- pc.res$metadata[,input$shape, drop = TRUE]
        if (is.factor(pl.shapes)) {
          pl.shapes <- droplevels(pl.shapes)
        }
      }

      if (!is.null(annot.by)) {
        hov <- list()
        for (i in seq_along(annot.by)) {
          hov[[i]] <- paste0("</br><b>",annot.by[i], ":</b> ", pc.res$metadata[[annot.by[i]]])
        }

        hov.text <- do.call(paste0, hov)
      }

      # Check if 2D is wanted.
      if (input$twod) {
        plot_ly(pc.res$rotated, x = as.formula(paste0("~",input$dim1)),
                y = as.formula(paste0("~",input$dim2)),
                type = "scatter",
                mode = "markers",
                marker = list(size = 15),
                color = pl.cols,
                colors = pl.col,
                symbol = pl.shapes,
                symbols = c("circle", "square", "diamond", "cross", "diamond-open", "circle-open", "square-open", "x"),
                text = hov.text,
                width = 1000, height = 710, hoverinfo = "text") %>%
          layout(xaxis = list(showgrid = FALSE, showline = TRUE, mirror = TRUE, zeroline = FALSE),
                 yaxis = list(showgrid = FALSE, showline = TRUE, mirror = TRUE, zeroline = FALSE))

      } else {

        # Generate plot.
        plot_ly(pc.res$rotated, x = as.formula(paste0("~",input$dim1)),
                y = as.formula(paste0("~",input$dim2)),
                z = as.formula(paste0("~",input$dim3)),
                type = "scatter3d",
                mode = "markers",
                color = pl.cols,
                colors = pl.col,
                symbol = pl.shapes,
                symbols = c("circle", "square", "diamond", "cross", "diamond-open", "circle-open", "square-open", "x"),
                text = hov.text,
                width = 1000, height = 710, hoverinfo = "text")
      }
    })
  }
  shinyApp(ui, server, options = list(height = height))
}
