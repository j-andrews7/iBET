#' Render general outputs for the Metadata tab
#'
#' Create rendering expressions for the Metadata outputs.
#'
#' @param input The Shiny input object from the server function.
#' @param output The Shiny output object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return A \linkS4class{NULL} is invisibly returned
#' and rendering expressions for Metadata tab features are added to \code{output}.
#'
#' @author Jared Andrews
#'
#' @importFrom shiny isolate selectInput tagList renderUI
#' @importFrom plotly renderPlotly %>%
#' @importFrom DT renderDT datatable formatStyle
#' @rdname INTERNAL_create_fgsea_metadata_outputs
.create_fgsea_metadata_outputs <- function(input, output, robjects) {
    # nocov start
    # output$fgsea.se.options <- renderUI({
    #     req(robjects$se)
    #     df <- robjects$se
    #     # Get only numeric variables.
    #     choices <- names(df)[vapply(df, is.numeric, logical(1))]

    #     # This prevents the selected value from resetting each time plots are re-rendered.
    #     # TODO: This seems gross, probably make a new function to handle UI output so it's not re-rendered
    #     # on plot updates.
    #     if (is.null(input$gene.esterm)) {
    #         esterm.sel <- ifelse("LFC" %in% choices, "LFC",
    #             ifelse("beta" %in% choices, "beta", choices[1])
    #         )
    #     } else {
    #         esterm.sel <- isolate(input$gene.esterm)
    #     }

    #     if (is.null(input$gene.sigterm)) {
    #         sigterm.sel <- ifelse("fdr" %in% choices, "fdr",
    #             ifelse("pval" %in% choices, "pval", choices[1])
    #         )
    #     } else {
    #         sigterm.sel <- isolate(input$gene.sigterm)
    #     }

    #     tagList(
    #         column(
    #             6,
    #             selectInput("gene.esterm", "Effect size term:",
    #                 choices = choices,
    #                 selected = esterm.sel
    #             )
    #         ),
    #         column(
    #             6,
    #             selectInput("gene.sigterm", "Significance term:",
    #                 choices = choices,
    #                 selected = sigterm.sel
    #             )
    #         )
    #     )
    # })
    # nocov end

    # nocov start
    output$fgsea.metadata <- renderDT(server = FALSE, {
        req(robjects$se)

        df <- colData(robjects$se)

        datatable(df,
            rownames = FALSE,
            filter = "top",
            extensions = c("Buttons"),
            caption = paste0(input$gene.sel1, " Gene Summary"),
            options = list(
                search = list(regex = TRUE),
                pageLength = 10,
                dom = "Blfrtip",
                buttons = c("copy", "csv", "excel", "pdf", "print")
            )
        ) %>% formatStyle(0, target = "row", lineHeight = "50%")
    })

    # nocov end
    invisible(NULL)
}