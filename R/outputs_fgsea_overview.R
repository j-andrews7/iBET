#' Render general outputs for the fgsea Overview tab
#'
#' Create rendering expressions for the fgsea Overview tab outputs.
#'
#' @param input The Shiny input object from the server function.
#' @param output The Shiny output object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return A \linkS4class{NULL} is invisibly returned
#'   and rendering expressions for the fgsea Overview tab features
#'   are added to \code{output}.
#'
#' @author Jared Andrews
#'
#' @importFrom shiny renderUI renderPlot tagList column selectInput isolate fluidRow
#' @importFrom plotly renderPlotly ggplotly layout config plot_ly toWebGL add_segments add_annotations %>%
#' @importFrom ggplot2 scale_x_discrete guide_axis
#' @importFrom MAGeCKFlute MapRatesView
#' @importFrom shinyWidgets updatePickerInput
#' @importFrom shinyjs js
#' @importFrom dittoSeq dittoColors
#' @importFrom grid grid.newpage grid.text
#' @importFrom DT renderDT datatable formatStyle
#' @importFrom stats cor as.formula
#' @rdname INTERNAL_create_qc_output
.create_fgsea_overview_outputs <- function(input, output, robjects) {
    # nocov start
    output$pca.comps <- renderUI({
        req(robjects$pc)
        pcs <- robjects$pc

        tagList(
            fluidRow(
                column(4, selectInput("dim1", "Dim1:", choices = pcs$components, selected = "PC1")),
                column(4, selectInput("dim2", "Dim2:", choices = pcs$components, selected = "PC2")),
                if (length(pcs$components) > 2) {
                    column(4, selectInput("dim3", "Dim3:", choices = pcs$components, selected = "PC3"))
                }
            )
        )
    })
    # nocov end

    # nocov start
    output$ov.bar <- renderPlotly({
        req(input$gsea.set)

        input$ov.bar.update

        pca.res <- robjects$pc

        colorb <- NULL
        shapeb <- NULL

        if (isolate(input$bip.color) != "") {
            colorb <- isolate(input$bip.color)
        }

        if (isolate(input$bip.shape) != "") {
            shapeb <- isolate(input$bip.shape)
        }

        if ((isolate(input$bip.twod) | length(pca.res$components) < 3)) {
            dizzy <- NULL
        } else {
            dizzy <- isolate(input$dim3)
        }

        if (any(unlist(lapply(list(pca.res, isolate(input$dim1), isolate(input$dim2)), is.null)))) {
            fig <- .empty_plot("No counts provided.\nNo PCA performed.", plotly = TRUE)
        } else {
            fig <- plot_pca_biplot(
                pca.res = pca.res,
                dim.x = isolate(input$dim1),
                dim.y = isolate(input$dim2),
                dim.z = dizzy,
                color.by = colorb,
                shape.by = shapeb,
                hover.info = "Label",
                show.loadings = isolate(input$bip.loadings),
                n.loadings = isolate(input$bip.n.loadings),
                pt.size = isolate(input$pca.pt.size)
            )
        }

        robjects$plot.ov.bar <- fig

        fig
    })
    # nocov end

    # nocov start
    output$ov.swoosh <- renderPlotly({
        if (!is.null(robjects$count.summary)) {
            fig <- plot_bar(robjects$count.summary)
        } else {
            fig <- .empty_plot("No summary provided.\nNo gini values available.", plotly = TRUE)
        }
        robjects$plot.ov.swoosh <- fig
        fig
    })
    # nocov end

    # nocov start
    output$ov.le <- renderPlotly({
        n.counts <- robjects$norm.counts
        n.counts.log <- as.matrix(log2(n.counts[, c(-1, -2)] + 1))
        colnames(n.counts.log) <- colnames(n.counts)[c(-1, -2)]

        # TODO: Add input to control gridlines.
        if (!is.null(n.counts.log) & nrow(n.counts.log) > 0) {
            fig <- plot_hist(n.counts.log,
                title = "Distribution of read counts",
                xlab = "log2(counts + 1)", ylab = "Frequency", show.grid = FALSE
            )
        } else {
            fig <- .empty_plot("No counts provided.\nCannot make distribution.", plotly = TRUE)
        }

        robjects$plot.ove.le <- fig
        fig
    })
    # nocov end

    # nocov start
    output$fgsea.summary <- renderDT(server = FALSE, {
        datatable(robjects$count.summary,
            rownames = FALSE,
            filter = "top",
            extensions = c("Buttons", "Scroller"),
            options = list(
                search = list(regex = TRUE),
                lengthMenu = list(c(10, 25, 50, -1), c("10", "25", "50", "all")),
                dom = "Blfrtip",
                buttons = c("copy", "csv", "excel", "pdf", "print"),
                scrollX = TRUE,
                deferRender = TRUE,
                scrollY = 600,
                scroller = TRUE
            )
        ) %>% formatStyle(0, target = "row", lineHeight = "80%")
    })
    # nocov end
    invisible(NULL)
}