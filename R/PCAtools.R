#' Create an interactive Shiny app for PCAtools analysis and results exploration
#'
#' @details Features with no variation will be removed prior to \code{\link[PCAtools]{pca}} being run.
#'
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @import DT
#' @import PCAtools
#' @importFrom plotly ggplotly plotlyOutput renderPlotly toWebGL layout plot_ly add_segments add_annotations
#' @import ggplot2
#' @import shinydashboard
#' @import dashboardthemes
#' @importFrom dittoSeq dittoColors
#' @importFrom grid grid.newpage grid.text
#' @importFrom shinyWidgets prettyCheckbox
#' @importFrom shinycustomloader withLoader
#' @importFrom shinyjqui jqui_resizable
#' @importFrom matrixStats rowVars
#' @importFrom stats as.formula p.adjust p.adjust.methods dist
#' @importFrom shinyBS bsCollapse bsCollapsePanel
#' @importFrom ComplexHeatmap pheatmap Heatmap HeatmapAnnotation
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
#' @seealso
#' \code{\link[PCAtools]{pca}}, \code{\link[PCAtools]{screeplot}},
#' \code{\link[PCAtools]{biplot}}, \code{\link[PCAtools]{pairsplot}}, \code{\link[PCAtools]{eigencorplot}}.
#'
#' @author Jared Andrews
#' @export
shinyPCAtools <- function(mat, metadata, removeVar = 0.3, scale = FALSE,
                          center = TRUE, color.by = NULL, shape.by = NULL,
                          annot.by = NULL, height = 850) {

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
            margin-bottom: 3px;
            padding-bottom: 2px !important;
            padding-top: 2px !important;
            font-size: 10px;
            line-height: 1.1;
          }
          .well {
            padding: 5px;
            margin-bottom: 10px;
          }
          .form-control, .selectize-input {
            padding-bottom: 2px !important;
            padding-top: 2px !important;
            font-size: 10px;
            height: 24px;
            min-height: 24px;
            line-height: 1.1;
          }
          .control-label {
            font-size: 10px;
            margin-bottom: 2px;
          }
          .panel-heading {
            padding: 5px 10px;
          }
          .selectize-control {
            margin-bottom: 0px;
          }
          body {
            line-height: 1.1;
          }
        "))
      ),
      shinyDashboardThemes(
        theme = "onenote"
      ),
      sidebarLayout(
        sidebarPanel(width = 3,
          bsCollapse(open = "biplot.settings",
            bsCollapsePanel(
              title = span(icon("plus"), "PCA Settings"), value = "pca.settings", style = "info",
              conditionalPanel(
                condition = "input['keep.top.n'] == false",
                numericInput("var.remove", "Remove this proportion of features ranked by variance:",
                             min = 0, max = 1, step = 0.01, value = removeVar)
              ),
              conditionalPanel(
                condition = "input['keep.top.n'] == true",
                numericInput("var.n.keep", "Number of features to retain by variance:",
                            min = 2, max = Inf, step = 1, value = 500)
              ),
              fluidRow(
                column(6,
                       prettyCheckbox("center", strong("Center data"), center, bigger = FALSE,
                                      animation = "smooth", status = "success",
                                      icon = icon("check"), width = "100%"),
                       prettyCheckbox("keep.top.n", strong("Limit by top N features"), FALSE, bigger = FALSE,
                                      animation = "smooth", status = "success",
                                      icon = icon("check"), width = "100%")
                ),
                column(6,
                       prettyCheckbox("scale", strong("Scale data"), scale, bigger = FALSE,
                                      animation = "smooth", status = "success",
                                      icon = icon("check"), width = "100%")
                )
              ),
              prettyCheckbox("meta.filt", strong("Filter via metadata table"), TRUE, bigger = FALSE,
                             animation = "smooth", status = "success",
                             icon = icon("check"), width = "100%")
            ),
            bsCollapsePanel(title = span(icon("plus"), "Biplot Settings"), value = "biplot.settings", style = "info",
              uiOutput("pca.comps"),
              fluidRow(
                column(6, selectInput("bip.color", "Color by:", choices = c("", colnames(metadata)),
                                      selected = color.by)),
                column(6, selectInput("bip.shape", "Shape by:", choices = c("", colnames(metadata)),
                                      selected = shape.by))
              ),
              fluidRow(
                column(6,
                       prettyCheckbox("bip.twod", strong("Limit to 2D"), FALSE, bigger = FALSE,
                         animation = "smooth", status = "success",
                         icon = icon("check"), width = "100%")
                ),
                column(6,
                       prettyCheckbox("bip.loadings", strong("Plot Loadings"), FALSE, bigger = FALSE,
                         animation = "smooth", status = "success",
                         icon = icon("check"), width = "100%")
                )
              ),
              numericInput("bip.n.loadings", "Loadings:",
                           min = 0, max = 100, step = 1, value = 5)
            ),
            bsCollapsePanel(title = span(icon("plus"), "Screeplot Settings"), value = "scree.settings", style = "info",
              fluidRow(
                column(6,
                       numericInput("scree.components", "Components:",
                                    min = 1, max = 50, step = 1, value = 30),
                       colourInput("scree.bar.col", "Bar color:", value = "#1E90FF"),
                       colourInput("scree.sumline.col", "Sum line color:", value = "#EE0000"),
                       numericInput("scree.sumline.cex", "Sum line size:", value = 1.5, min = 0.01, step = 0.1),
                       colourInput("scree.sumpts.col", "Sum points color:", value = "#EE0000"),
                       numericInput("scree.sumpts.cex", "Sum points size:", value = 2, min = 0.01, step = 0.1)
                ),
                column(6,
                       textInput("scree.main", "Main title:", value = "", placeholder = "Enter text"),
                       numericInput("scree.hline.val", "Hline value:", value = 80, min = 0.01, step = 0.5),
                       prettyCheckbox("scree.horns", strong("Plot Horn's parallel"), FALSE, bigger = FALSE,
                                      animation = "smooth", status = "success",
                                      icon = icon("check"), width = "100%"),
                       prettyCheckbox("scree.elbow", strong("Plot elbow point"), FALSE, bigger = FALSE,
                                      animation = "smooth", status = "success",
                                      icon = icon("check"), width = "100%"),
                       prettyCheckbox("scree.sumline", strong("Draw sum line"), TRUE, bigger = FALSE,
                                      animation = "smooth", status = "success",
                                      icon = icon("check"), width = "100%"),
                       prettyCheckbox("scree.sumpts", strong("Draw sum points"), TRUE, bigger = FALSE,
                                      animation = "smooth", status = "success",
                                      icon = icon("check"), width = "100%"),
                       prettyCheckbox("scree.grid.maj", strong("Draw gridlines"), TRUE, bigger = FALSE,
                                      animation = "smooth", status = "success",
                                      icon = icon("check"), width = "100%"),
                       prettyCheckbox("scree.hline", strong("Draw horizontal line"), FALSE, bigger = FALSE,
                                      animation = "smooth", status = "success",
                                      icon = icon("check"), width = "100%")
                )
              )
            ),
            bsCollapsePanel(title = span(icon("plus"), "Eigencorplot Settings"), value = "eigen.settings", style = "info",
              fluidRow(
                column(6,
                  numericInput("eig.components", "Components:",
                               min = 1, max = 50, step = 1, value = 10),
                  selectInput("eig.corfun", "Correlation:", choices = c("pearson", "spearman", "kendall"),
                              selected = "pearson"),
                  colourInput("eig.max.col", "Max color:", value = "#8B0000"),
                  colourInput("eig.mid.col", "Midpoint color:", value = "#FFFFFF"),
                  colourInput("eig.min.col", "Min color:", value = "#00008B"),
                  textInput("eig.main", "Main title:", value = "", placeholder = "Enter text"),
                  numericInput("eig.main.cex", "Main title size:", value = 2, min = 0.01, step = 0.1),
                  selectInput("eig.main.style", "Main title style:",
                              choices = list("plain", "bold", "italic"),
                              selected = "bold"),
                  numericInput("eig.corr.cex", "Corr values size:", value = 1, min = 0.01, step = 0.1),
                  selectInput("eig.corr.style", "Corr values style:",
                              choices = list("plain", "bold", "italic"),
                              selected = "plain"),
                  colourInput("eig.corr.col", "Corr values color:", value = "#000000")
                ),
                column(6,
                  uiOutput("eigen.vars"),
                  prettyCheckbox("eig.rsquare", strong("Plot R Square"), FALSE, bigger = FALSE,
                                 animation = "smooth", status = "success",
                                 icon = icon("check"), width = "100%"),
                  selectInput("eig.cormulttest", "Multiple Test Correction:",
                              choices = p.adjust.methods, selected = "none"),
                  textInput("eig.x", "X title:", value = "", placeholder = "Enter text"),
                  numericInput("eig.x.cex", "X title size:", value = 1, min = 0.01, step = 0.1),
                  selectInput("eig.x.style", "X title style:",
                              choices = list("plain", "bold", "italic")),
                  textInput("eig.y", "Y title:", value = "", placeholder = "Enter text"),
                  numericInput("eig.y.cex", "Y title size:", value = 1, min = 0.01, step = 0.1),
                  selectInput("eig.y.style", "Y title style:",
                              choices = list("plain", "bold", "italic")),
                  selectInput("eig.posKey", "Key position:",
                              choices = list("right", "left", "top", "bottom"))
                )
              )
            ),
            bsCollapsePanel(title = span(icon("plus"), "Distance Matrix Settings"), value = "dist.settings", style = "info",
              fluidRow(
                column(6,
                  colourInput("dist.min.col", "Min color:", value = "darkblue"),
                  colourInput("dist.max.col", "Max color:", value = "#FFFFFF"),
                  numericInput("dist.row.cex", "Row font size:", value = 10, min = 1, step = 0.25)
                ),
                column(6,
                  selectInput("dist.method", "Method:",
                              choices = c("euclidean", "maximum", "manhattan", "canberra", "minkowski"),
                              selected = "euclidean"),
                  prettyCheckbox("dist.rownames", label = "Show row names", value = TRUE,
                                 animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                  prettyCheckbox("dist.colnames", label = "Show column names", value = FALSE,
                                 animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                  numericInput("dist.col.cex", "Col font size:", value = 10, min = 1, step = 0.25)
                )
              ),
              fluidRow(
                uiOutput("dist.anno.opts")
              )
            )
          ),
          div(actionButton("update", "Update Plots"), align = "center")
        ),
        mainPanel(width = 9,
          tabsetPanel(
            tabPanel("biplot", div(jqui_resizable(plotlyOutput("biplot", height = "700px", width = "1000px")), align = "center")),
            tabPanel("screeplot", div(jqui_resizable(plotlyOutput("screeplot")), align = "center")),
            tabPanel("eigencorplot", div(jqui_resizable(plotOutput("eigencorplot")), align = "center")),
            tabPanel("Distance Matrix", div(jqui_resizable(plotOutput("distmatrix")), align = "center")),
            tabPanel("Metadata (Filtering)", div(br(), DTOutput("metadata"), style = "font-size:80%"))
          )
        )
      )
    )
  )

  server <- function(input, output, session) {

    matty <- reactive({
      matt <- mat

      if (!is.null(input$metadata_rows_all) & input$meta.filt) {
        matt <- matt[,input$metadata_rows_all]
      }

      # Remove features with no variance.
      matt <- matt[rowVars(matt) > 0,]

      # If necessary, limit to top N features by variance.
      if (input$keep.top.n) {
        matt <- matt[order(rowVars(matt), decreasing = TRUE),]

        if (input$var.n.keep < nrow(matt)) {
          matt <- matt[1:input$var.n.keep,]
        }
      }

      matt
    })

    pc <- reactive({
      req(input$var.remove)
      meta <- metadata
      mat <- matty()

      if (!is.null(input$metadata_rows_all) & input$meta.filt) {
        meta <- metadata[input$metadata_rows_all,]
      }

      # If input to use top N features instead rather than percent-based feature removal, account for that
      if (input$keep.top.n) {
        var.remove <- 0
        mat <- mat[order(rowVars(mat), decreasing = TRUE),]
        mat <- mat[1:input$keep.top.n,]
      } else {
        var.remove <- input$var.remove
      }

      pca(mat, metadata = meta, removeVar = var.remove, scale = input$scale, center = input$center)
    })

    nonnum_vars <- reactive({
      req(pc)

      pcs <- pc()
      pcs$metadata[,!unlist(lapply(pcs$metadata, is.numeric))]
    })

    output$dist.anno.opts <- renderUI({
      req(nonnum_vars)
      local({
        nonnum_var <- nonnum_vars()

        tagList(
          column(6,
                 selectInput("dist.top.anno", "Top annotation:", choices = c("", names(nonnum_var)), multiple = TRUE),
                 selectInput("dist.bot.anno", "Bottom annotation:", choices = c("", names(nonnum_var)), multiple = TRUE)
          ),
          column(6,
                 selectInput("dist.right.anno", "Right annotation:", choices = c("", names(nonnum_var)), multiple = TRUE),
                 selectInput("dist.left.anno", "Left annotation:", choices = c("", names(nonnum_var)), multiple = TRUE)
          )
        )
      })
    })

    # Populate UI with all PCs.
    # To do: Write check for only 2 PCs.
    output$pca.comps <- renderUI({
      req(pc)
      local({
        pcs <- pc()

        tagList(
          fluidRow(
            column(4, selectInput("dim1", "Dim1:", choices = pcs$components, selected = "PC1")),
            column(4, selectInput("dim2", "Dim2:", choices = pcs$components, selected = "PC2")),
            column(4, selectInput("dim3", "Dim3:", choices = pcs$components, selected = "PC3"))
          )
        )
      })
    })

    # Populate eigencorplot UI with only the numeric metadata variables as choices.
    output$eigen.vars <- renderUI({
      req(pc)
      local({
        pcs <- pc()

        mets <- pcs$metadata[,unlist(lapply(pcs$metadata, is.numeric))]

        tagList(
          pickerInput("eig.vars", "Variables:", choices = c("", names(mets)),
                      multiple = TRUE, options = list(`live-search` = TRUE,
                                                      `actions-box` = TRUE))
        )
      })
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

    output$biplot <- renderPlotly({
      req(pc, input$dim1, input$dim2, input$dim3)
      input$update

      pc.res <- isolate(pc())

      pl.cols <- NULL
      pl.shapes <- NULL
      pl.col <- "black"
      hov.text <- NULL

      # Get marker aesthetics mappings.
      # Drop unused factor levels if possible.
      if (isolate(input$bip.color) != "") {
        pl.cols <- pc.res$metadata[,isolate(input$bip.color), drop = TRUE]
        if (is.factor(pl.cols)) {
          pl.cols <- droplevels(pl.cols)
        }
        pl.col <- dittoColors()[seq_along(unique(pc.res$metadata[,isolate(input$bip.color), drop = TRUE]))]
      }

      if (isolate(input$bip.shape) != "") {
        pl.shapes <- pc.res$metadata[,isolate(input$bip.shape), drop = TRUE]
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
      if (isolate(input$bip.twod)) {
        fig <- plot_ly(pc.res$rotated, x = as.formula(paste0("~", isolate(input$dim1))),
                y = as.formula(paste0("~", isolate(input$dim2))),
                type = "scatter",
                mode = "markers",
                marker = list(size = 15),
                color = pl.cols,
                colors = pl.col,
                symbol = pl.shapes,
                symbols = c("circle", "square", "diamond", "cross",
                            "diamond-open", "circle-open", "square-open", "x"),
                text = hov.text,
                hoverinfo = "text") %>%
          layout(xaxis = list(showgrid = FALSE, showline = TRUE, mirror = TRUE, zeroline = FALSE,
                              title = paste0(isolate(input$dim1),
                                             " (", format(round(pc.res$variance[isolate(input$dim1)], 2), nsmall = 2),"%)")),
                 yaxis = list(showgrid = FALSE, showline = TRUE, mirror = TRUE, zeroline = FALSE,
                              title = paste0(isolate(input$dim2),
                                             " (", format(round(pc.res$variance[isolate(input$dim2)], 2), nsmall = 2),"%)")))

        fig <- fig %>% toWebGL()

        # Plot loadings.
        if (isolate(input$bip.loadings)) {
          lengthLoadingsArrowsFactor <- 1.5

          # Get number of loadings to display.
          xidx <- order(abs(pc.res$loadings[,isolate(input$dim1)]), decreasing = TRUE)
          yidx <- order(abs(pc.res$loadings[,isolate(input$dim2)]), decreasing = TRUE)
          vars <- unique(c(
            rownames(pc.res$loadings)[xidx][seq_len(isolate(input$bip.n.loadings))],
            rownames(pc.res$loadings)[yidx][seq_len(isolate(input$bip.n.loadings))]))

          # get scaling parameter to match between variable loadings and rotated loadings
          # This is cribbed almost verbatim from PCAtools code.
          r <- min(
            (max(pc.res$rotated[,isolate(input$dim1)]) - min(pc.res$rotated[,isolate(input$dim1)]) /
               (max(pc.res$loadings[,isolate(input$dim1)]) - min(pc.res$loadings[,isolate(input$dim1)]))),
            (max(pc.res$rotated[,isolate(input$dim2)]) - min(pc.res$rotated[,isolate(input$dim2)]) /
               (max(pc.res$loadings[,isolate(input$dim2)]) - min(pc.res$loadings[,isolate(input$dim2)]))))

          fig <- fig %>%
            add_segments(x = 0, xend = pc.res$loadings[vars,isolate(input$dim1)] * r * lengthLoadingsArrowsFactor,
                         y = 0, yend = pc.res$loadings[vars,isolate(input$dim2)] * r * lengthLoadingsArrowsFactor,
                         line = list(color = 'black'), inherit = FALSE, showlegend = FALSE, hoverinfo = "text") %>%
            add_annotations(x = pc.res$loadings[vars,isolate(input$dim1)] * r * lengthLoadingsArrowsFactor,
                            y = pc.res$loadings[vars,isolate(input$dim2)] * r * lengthLoadingsArrowsFactor,
                            ax = 0, ay = 0, text = vars, xanchor = 'center', yanchor= 'bottom')
        }
      } else {

        # Generate plot.
        fig <- plot_ly(pc.res$rotated, x = as.formula(paste0("~", isolate(input$dim1))),
                y = as.formula(paste0("~", isolate(input$dim2))),
                z = as.formula(paste0("~", isolate(input$dim3))),
                type = "scatter3d",
                mode = "markers",
                color = pl.cols,
                colors = pl.col,
                symbol = pl.shapes,
                symbols = c("circle", "square", "diamond", "cross", "diamond-open",
                            "circle-open", "square-open", "x"),
                text = hov.text,
                hoverinfo = "text") %>%
          layout(scene = list(
            xaxis = list(title = paste0(isolate(input$dim1), " (",
                                        format(round(pc.res$variance[isolate(input$dim1)], 2), nsmall = 2),"%)")),
            yaxis = list(title = paste0(isolate(input$dim2), " (",
                                        format(round(pc.res$variance[isolate(input$dim2)], 2), nsmall = 2),"%)")),
            zaxis = list(title = paste0(isolate(input$dim3), " (",
                                        format(round(pc.res$variance[isolate(input$dim3)], 2), nsmall = 2),"%)")),
            camera = list(eye = list(x=1.5, y = 1.8, z = 0.4))))
      }
      fig <- fig %>%
        config(edits = list(annotationPosition = TRUE,
                            annotationTail = FALSE),
               toImageButtonOptions = list(format = "svg"),
               displaylogo = FALSE,
               plotGlPixelRatio = 7)

      fig
    })

    output$screeplot <- renderPlotly({
      req(pc)
      input$update

      pc.res <- isolate(pc())

      # Limit the components to those that exist.
      comps <- ifelse(isolate(input$scree.components) > length(pc.res$components),
                      length(pc.res$components), isolate(input$scree.components))

      horn <- NULL
      if (isolate(input$scree.horns)) {
        horn <- parallelPCA(matty())
        horn <- horn$n
      }

      elbow <- NULL
      if (isolate(input$scree.elbow)) {
        elbow <- findElbowPoint(pc.res$variance)
      }

      # Check/get hline value.
      hline <- NULL
      if (isolate(input$scree.hline)) {
        hline <- isolate(input$scree.hline.val)
      }

      # Make plot.
      gg <- screeplot(pc.res,
                components = getComponents(pc.res, seq_len(comps)),
                colBar = isolate(input$scree.bar.col),
                colCumulativeSumLine = isolate(input$scree.sumline.col),
                sizeCumulativeSumLine = isolate(input$scree.sumline.cex),
                colCumulativeSumPoints = isolate(input$scree.sumpts.col),
                sizeCumulativeSumPoints = isolate(input$scree.sumpts.cex),
                title = isolate(input$scree.main),
                drawCumulativeSumLine = isolate(input$scree.sumline),
                drawCumulativeSumPoints = isolate(input$scree.sumpts),
                gridlines.major = isolate(input$scree.grid.maj),
                hline = hline,
                vline = c(horn, elbow),
                xlabAngle = 45)

      fig <- ggplotly(gg, tooltip = c("x", "y")) %>%
        config(edits = list(annotationPosition = TRUE,
                            annotationTail = FALSE),
               toImageButtonOptions = list(format = "svg"),
               displaylogo = FALSE,
               plotGlPixelRatio = 7)

      # Add vline annotations if plotted.
      if (!is.null(horn)) {
        fig <- fig %>% add_annotations(x = c(horn),
                         y = c(50),
                         text = c("Horn's"),
                         showarrow = FALSE,
                         font = list(size = 16))
      }

      if (!is.null(elbow)) {
        fig <- fig %>% add_annotations(x = c(elbow),
                                       y = c(50),
                                       text = c("Elbow point"),
                                       showarrow = FALSE,
                                       font = list(size = 16))
      }

      fig
    })

    output$eigencorplot <- renderPlot({
      req(pc)
      input$update

      pc.res <- isolate(pc())

      # Limit the components to those that exist.
      comps <- ifelse(isolate(input$eig.components) > length(pc.res$components),
                      length(pc.res$components), isolate(input$eig.components))

      if (!is.null(isolate(input$eig.vars)) & length(isolate(input$eig.vars)) > 1) {
        eigencorplot(pcaobj = pc.res,
                     components = getComponents(pc.res, seq_len(comps)),
                     metavars = isolate(input$eig.vars),
                     col = c(isolate(input$eig.min.col), isolate(input$eig.mid.col), isolate(input$eig.max.col)),
                     plotRsquared = isolate(input$eig.rsquare),
                     corMultipleTestCorrection = isolate(input$eig.cormulttest),
                     corFUN = isolate(input$eig.corfun),
                     main = isolate(input$eig.main),
                     cexMain = isolate(input$eig.main.cex),
                     fontMain = isolate(input$eig.main.style),
                     titleX = isolate(input$eig.x),
                     cexTitleX = isolate(input$eig.x.cex),
                     fontTitleX = isolate(input$eig.x.style),
                     titleY = isolate(input$eig.y),
                     cexTitleY = isolate(input$eig.y.cex),
                     fontTitleY = isolate(input$eig.y.style),
                     rotTitleY = 90,
                     rotLabX = 45,
                     cexCorval = isolate(input$eig.corr.cex),
                     fontCorval = isolate(input$eig.corr.style),
                     colCorval = isolate(input$eig.corr.col),
                     posColKey = isolate(input$eig.posKey))
      } else {
        grid.newpage()
        grid.text("Select at least two numeric metadata values.")
      }
    })

    output$distmatrix <- renderPlot({
      req(matty)
      input$update

      ds.colors <- dittoColors()

      # Get metadata for samples remaining.
      meta <- metadata
      if (!is.null(isolate(input$metadata_rows_all))) {
        meta <- metadata[isolate(input$metadata_rows_all),]
      }

      left.anno <- NULL
      right.anno <- NULL
      top.anno <- NULL
      bot.anno <- NULL

      if (!is.null(isolate(input$dist.left.anno))) {
        left.anno <- .create_anno(input$dist.left.anno, meta, ds.colors, anno_type = "row", side = "top")
      }

      if (!is.null(isolate(input$dist.right.anno))) {
        right.anno <- .create_anno(input$dist.right.anno, meta, ds.colors, anno_type = "row", side = "bottom")
      }

      if (!is.null(isolate(input$dist.top.anno))) {
        top.anno <- .create_anno(input$dist.top.anno, meta, ds.colors, anno_type = "column", side = "right")
      }

      if (!is.null(isolate(input$dist.bot.anno))) {
        bot.anno <- .create_anno(input$dist.bot.anno, meta, ds.colors, anno_type = "column", side = "left")
      }

      dists <- dist(t(matty()), method = isolate(input$dist.method))
      sampleDistMatrix <- as.matrix(dists)

      colors <- c(isolate(input$dist.min.col), isolate(input$dist.max.col))

      ComplexHeatmap::pheatmap(sampleDistMatrix,
               clustering_distance_rows = dists,
               clustering_distance_cols = dists,
               col = colors,
               row_km = isolate(input$dist.row.km),
               column_km = isolate(input$dist.col.km),
               fontsize_row = isolate(input$dist.row.cex),
               fontsize_col = isolate(input$dist.col.cex),
               show_colnames = isolate(input$dist.colnames),
               show_rownames = isolate(input$dist.rownames),
               left_annotation = left.anno,
               right_annotation = right.anno,
               top_annotation = top.anno,
               bottom_annotation = bot.anno,
               heatmap_legend_param = list(title = paste0(isolate(input$dist.method), " Distance")))
    })

    # Initialize plots by simulating button click once.
    o <- observe({
      req(pc, input$dim1, input$dim2, input$dim3)
      shinyjs::click("update")
      o$destroy
    })
  }
  shinyApp(ui, server, options = list(height = height))
}
